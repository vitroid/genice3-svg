from logging import getLogger
from math import atan2, exp, pi
import colorsys


import svgwrite as sw
import numpy as np


def polygon_path(vs, **kwargs):
    p = []
    p.append(["M", vs[-1][0], vs[-1][1]])
    for v in vs:
        p.append(["L", v[0], v[1]])
    p.append(["Z"])
    return sw.path.Path(d=p, **kwargs)


def polygon(svg, com, d, **options):
    """
    draw a polygon
    """
    group = svg.add(svg.g(id="Polygon"))
    path = polygon_path(com + d, **options)
    group.add(path)


sun = np.array([-1.0, -1.0, 2.0])
sun /= np.linalg.norm(sun)


def Normal(vs):
    """
    Normal vector (not normalized)
    """
    n = np.zeros(3)
    for i in range(vs.shape[0]):
        n += np.cross(vs[i - 1], vs[i])
    return n


def cylinder_path(R, ratio, L, **kwargs):
    # horizontal, start from origin
    magic = 0.552284749831
    x1 = R * ratio
    x2 = x1 * magic
    y1 = R
    y2 = y1 * magic
    p = []
    p.append(["M", 0, -y1])
    p.append(["L", L, -y1])
    p.append(["C", L + x2, -y1, L + x1, -y2, L + x1, 0])
    p.append(["C", L + x1, y2, L + x2, y1, L, y1])
    p.append(["L", 0, y1])
    p.append(["C", -x2, y1, -x1, y2, -x1, 0])
    p.append(["C", -x1, -y2, -x2, -y1, 0, -y1])
    p.append(["Z"])
    return sw.path.Path(d=p, **kwargs)


def cylinder_new(svg, v1_, v2_, r, **options):
    """
    draw a 3D cylinder
    """
    logger = getLogger()
    group = svg.add(svg.g(id="Cylinder"))
    if v1_[2] > v2_[2]:
        v1, v2 = v2_, v1_
    else:
        v1, v2 = v1_, v2_
    dir = v2[:2] - v1[:2]
    angle = atan2(dir[1], dir[0])
    d = v2 - v1
    ratio = d[2] / np.linalg.norm(d)
    L = np.linalg.norm(dir)
    path = cylinder_path(r, ratio, L, **options)
    path.translate(v1[0], v1[1])
    path.rotate(angle * 180 / pi, center=(0, 0))
    group.add(path)
    u = sw.shapes.Ellipse(center=v2[:2], r=(ratio * r, r), **options)  # , fill=endfill)
    u.rotate(angle * 180 / pi, center=v2[:2])
    group.add(u)


def Render(
    prims,
    Rsphere,
    shadow=None,
    zoom=100,
    topleft=np.array([-1.0, -1.0]),
    size=(50.0, 50.0),
    bgcolor=None,
    encode=False,
):
    # encode is dummy parameter.
    logger = getLogger()
    svg = sw.Drawing(
        size=("{0}px".format(size[0] * zoom), "{0}px".format(size[1] * zoom))
    )
    if bgcolor is not None:
        svg.add(
            svg.rect(
                insert=(0, 0),
                size=("{0}px".format(size[0] * zoom), "{0}px".format(size[1] * zoom)),
                fill=bgcolor,
            )
        )
    TL0 = np.zeros(3)
    TL0[:2] = topleft
    shadows = []
    linedefaults = {
        "stroke_width": 2,
        "stroke": "#000",
        "stroke_linejoin": "round",
        "stroke_linecap": "round",
    }
    filldefaults = {
        "stroke_width": 1,
        "stroke": "#000",
        "fill": "#0ff",
        "stroke_linejoin": "round",
        "stroke_linecap": "round",
        "fill_opacity": 1.0,
    }
    shadowdefaults = {
        "stroke_width": 0,
        "fill": "#888",
        "fill_opacity": 0.08,
    }
    if shadow is not None:
        assert shadow[0] == "#" and len(shadow) in (5, 9)
        if len(shadow) == 5:
            shadowdefaults["fill"] = shadow[:4]
            shadowdefaults["fill_opacity"] = int(shadow[4], 16) / 15
        else:
            shadowdefaults["fill"] = shadow[:7]
            shadowdefaults["fill_opacity"] = int(shadow[7:], 16) / 255
        r = Rsphere
        Z = np.array([0, 0, 1.0])
        prims += [
            [prim[0] - Z * r * 1.4**j, prim[1] + "S", r * 1.4**j, {}]
            for j in range(1, 5)
            for prim in prims
            if prim[1] == "C"
        ]
    for prim in sorted(prims, key=lambda x: x[0][2]):
        if not (
            (-0.5 < prim[0][0] - topleft[0] < size[0] + 0.5)
            and (-0.5 < prim[0][1] - topleft[1] < size[1] + 0.5)
        ):
            continue
        if prim[1] == "L":
            if prim[4] == 0:
                options = {**linedefaults, **prim[5]}
                svg.add(
                    sw.shapes.Line(
                        start=(prim[2][:2] - topleft) * zoom,
                        end=(prim[3][:2] - topleft) * zoom,
                        **options,
                    )
                )
            else:
                options = {**filldefaults, **prim[5]}
                cylinder_new(
                    svg,
                    (prim[2] - TL0) * zoom,
                    (prim[3] - TL0) * zoom,
                    prim[4] * zoom,
                    **options,
                )
        elif prim[1] == "L2":
            if prim[3] == 0:
                options = {**linedefaults, **prim[4]}
                s = ((prim[0] + prim[2])[:2] - topleft) * zoom
                e = ((prim[0] - prim[2])[:2] - topleft) * zoom
                svg.add(sw.shapes.Line(start=s, end=e, **options))
            else:
                options = {**filldefaults, **prim[4]}
                cylinder_new(
                    svg,
                    (prim[0] + prim[2] - TL0) * zoom,
                    (prim[0] - prim[2] - TL0) * zoom,
                    prim[3] * zoom,
                    **options,
                )
        elif prim[1] == "P":
            options = prim[3]
            if "fillhs" in options:
                normal = Normal(prim[2])
                normal /= np.linalg.norm(normal)
                cosine = abs(np.dot(sun, normal))
                hue, sat = options["fillhs"]
                del options["fillhs"]
                bri = cosine * 0.5 + 0.5
                if sat < 0.2:
                    bri *= 0.9
                if cosine > 0.8:
                    sat *= 1 - (cosine - 0.8) * 3
                r, g, b = colorsys.hsv_to_rgb(hue / 360.0, sat, bri)
                rgb = "#{0:x}{1:x}{2:x}".format(
                    int(r * 15.9), int(g * 15.9), int(b * 15.9)
                )
                options["fill"] = rgb
            options = {**filldefaults, **options}
            polygon(svg, (prim[0] - TL0) * zoom, prim[2] * zoom, **options)
        elif prim[1] == "C":
            options = {**filldefaults, **prim[3]}
            Rsphere = prim[2]
            svg.add(
                sw.shapes.Circle(
                    center=(prim[0][:2] - topleft) * zoom, r=Rsphere * zoom, **options
                )
            )
        elif prim[1] == "CS":
            Rsphere = prim[2]
            options = {**shadowdefaults, **prim[3]}
            svg.add(
                sw.shapes.Circle(
                    center=(prim[0][:2] - topleft) * zoom, r=Rsphere * zoom, **options
                )
            )
    return svg.tostring()






