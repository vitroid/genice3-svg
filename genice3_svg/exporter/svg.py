# coding: utf-8
"""
GenIce exporter plugin to generate a SVG file.

Usage:
    % genice3 CS2 -r 3 3 3 -e svg[rotatex=30,shadow] > CS2.svg

Options:
    rotatex=30
    rotatey=30
    rotatez=30
    polygon        Draw polygons instead of a ball and stick model.
    arrows         Draw the hydrogen bonds with arrows.
    shadow=#8881   Draw shadows behind balls.
    bg=#f00        Specify the background color.
    O=0.06
    H=0            Size of the hydrogen atom (relative to that of oxygen)
    HB=0.4         Radius of HB relative to that of oxygem
    OH=0.5         Radius of OH colvalent bond relative to that of oxygem
    width=0        (Pixel)
    height=0       (Pixel)
    margin=0       (Pixel)
"""

# from genice2.molecules import serialize
# from genice2.decorators import timeit, banner
# import genice2.formats
from genice3_svg.render_svg import Render
from cycless.cycles import cycles_iter
import networkx as nx
import numpy as np
from collections import defaultdict
from logging import getLogger
from math import pi, cos, sin, radians
import re
from dataclasses import dataclass
from genice3.genice import GenIce3
import sys
from io import TextIOWrapper
from genice3.exporter import (
    parse_water_model_option,
)

desc = {
    "ref": {},
    "brief": "SVG (Standard Vector Graphics).",
    "usage": __doc__,
}


def Normal(vs):
    """
    Normal vector (not normalized)
    """
    n = np.zeros(3)
    for i in range(vs.shape[0]):
        n += np.cross(vs[i - 1], vs[i])
    return n


sun = np.array([-1.0, -1.0, 2.0])
sun /= np.linalg.norm(sun)


# set of hue and saturation
hue_sat = {
    3: (60.0, 0.8),
    4: (120, 0.8),  # yellow-green
    5: (180, 0.5),  # skyblue
    6: (240, 0.5),  # blue
    7: (300, 0.8),
    8: (350, 0.5),
}  # red-purple


def clip_cyl(v1, r1, v2, r2, rb):
    r1c = (r1**2 - rb**2) ** 0.5
    r2c = (r2**2 - rb**2) ** 0.5
    dv = v2 - v1
    Lv = np.linalg.norm(dv)
    if Lv < r1 + r2:
        return None
    newv1 = v1 + dv * r1c / Lv
    newv2 = v2 - dv * r2c / Lv
    c = (newv1 + newv2) / 2
    d = c - newv2
    return [c, "L2", d]


def draw_cell(prims, cellmat, origin=np.zeros(3)):
    for a in (0.0, 1.0):
        for b in (0.0, 1.0):
            v0 = np.array([0.0, a, b] + origin)
            v1 = np.array([1.0, a, b] + origin)
            mid = (v0 + v1) / 2
            prims.append(
                [
                    np.dot(mid, cellmat),
                    "L",
                    np.dot(v0, cellmat),
                    np.dot(v1, cellmat),
                    0,
                    {"stroke_width": 1, "stroke": "#888"},
                ]
            )
            v0 = np.array([b, 0.0, a] + origin)
            v1 = np.array([b, 1.0, a] + origin)
            mid = (v0 + v1) / 2
            prims.append(
                [
                    np.dot(mid, cellmat),
                    "L",
                    np.dot(v0, cellmat),
                    np.dot(v1, cellmat),
                    0,
                    {"stroke_width": 1, "stroke": "#888"},
                ]
            )
            v0 = np.array([a, b, 0.0] + origin)
            v1 = np.array([a, b, 1.0] + origin)
            mid = (v0 + v1) / 2
            prims.append(
                [
                    np.dot(mid, cellmat),
                    "L",
                    np.dot(v0, cellmat),
                    np.dot(v1, cellmat),
                    0,
                    {"stroke_width": 1, "stroke": "#888"},
                ]
            )
    corners = []
    for x in (np.zeros(3), cellmat[0]):
        for y in (np.zeros(3), cellmat[1]):
            for z in (np.zeros(3), cellmat[2]):
                corners.append(x + y + z + origin)
    corners = np.array(corners)
    return (
        np.min(corners[:, 0]),
        np.max(corners[:, 0]),
        np.min(corners[:, 1]),
        np.max(corners[:, 1]),
    )


from genice3.dictparser import parse_dict_options

type_map = dict(
    encode=bool,
    poly=bool,
    shadow=str,
    oxygen=float,
    HB=float,
    OH=float,
    hydrogen=float,
    arrows=bool,
    bgcolor=str,
    proj=np.array,
    width=int,
    height=int,
    margin=int,
    unprocessed=dict,
)

defaults = dict(
    encode=True,
    poly=False,
    shadow=None,
    oxygen=0.06,
    HB=0.4,
    OH=0.5,
    hydrogen=0,
    arrows=False,
    bgcolor=None,
    proj=np.array([[1.0, 0, 0], [0, 1, 0], [0, 0, 1]]),
    width=0,
    height=0,
    margin=0,
    unprocessed=dict(),
)

list_keys = dict(rotate=",")


def parse_options(**kwargs):
    options = parse_dict_options(
        kwargs, type_map=type_map, defaults=defaults, list_keys=list_keys
    )
    logger = getLogger()
    # encode = True  # valid for png.
    # poly = False
    # shadow = None
    # oxygen = 0.06  # absolute radius in nm
    # HB = 0.4  # radius relative to the oxygen
    # OH = 0.5  # radius relative to the oxygen
    # hydrogen = 0  # radius relative to the oxygen
    # arrows = False
    # bgcolor = None
    # proj = np.array([[1.0, 0, 0], [0, 1, 0], [0, 0, 1]])
    # width = 0
    # height = 0
    # margin = 0  # margin around the cell cube
    # unprocessed = dict()
    # for key, value in kwargs.items():
    #     logger.info(f"  Option with arguments: {key} := {value}")
    #     if key == "rotmat":
    #         value = re.search(r"\[([-0-9,.]+)\]", value).group(1)
    #         proj = np.array([float(x) for x in value.split(",")]).reshape(3, 3)
    #     elif key == "rotatex":
    #         logger.warning(
    #             "The rotatey option is deprecated. Use rotate option instead."
    #         )
    #         value = float(value) * pi / 180
    #         cosx = cos(value)
    #         sinx = sin(value)
    #         R = np.array([[cosx, 0, -sinx], [0, 1, 0], [sinx, 0, cosx]])
    #         proj = np.dot(proj, R)
    #     elif key == "rotatez":
    #         logger.warning(
    #             "The rotatez option is deprecated. Use rotate option instead."
    #         )
    #         value = float(value) * pi / 180
    #         cosx = cos(value)
    #         sinx = sin(value)
    #         R = np.array([[cosx, sinx, 0], [-sinx, cosx, 0], [0, 0, 1]])
    #         proj = np.dot(proj, R)
    #     elif key == "rotate":
    #         values = value.split(",")
    #         for value in values:
    #             axis = value[0]
    #             angle = radians(float(value[1:]))
    #             cosx = cos(angle)
    #             sinx = sin(angle)
    #             if axis in "xX":
    #                 R = np.array([[1, 0, 0], [0, cosx, sinx], [0, -sinx, cosx]])
    #             elif axis in "yY":
    #                 R = np.array([[cosx, 0, -sinx], [0, 1, 0], [sinx, 0, cosx]])
    #             elif axis in "zZ":
    #                 R = np.array([[cosx, sinx, 0], [-sinx, cosx, 0], [0, 0, 1]])
    #             else:
    #                 assert False, "  Wrong options."
    #             proj = np.dot(proj, R)
    #     elif key == "shadow":
    #         if value is True:
    #             shadow = "#8881"
    #         else:
    #             shadow = value
    #     elif key == "H":
    #         if value is True:
    #             hydrogen = 0.6
    #             HB = 0.2
    #         else:
    #             hydrogen = float(value)
    #     elif key == "HB":
    #         HB = float(value)
    #     elif key == "O":
    #         oxygen = float(value)
    #     elif key == "OH":
    #         if value is True:
    #             OH = 0.5
    #         else:
    #             OH = float(value)
    #     elif key == "bg":
    #         bgcolor = value
    #     elif key == "width":
    #         width = int(value)
    #     elif key == "height":
    #         height = int(value)
    #     elif key == "margin":
    #         margin = int(value)
    #     elif key == "encode":
    #         encode = bool(value)
    #     elif value is True:
    #         a = key
    #         logger.info("  Flags: {0}".format(a))
    #         if a == "polygon":
    #             poly = True
    #         elif a == "arrows":
    #             arrows = True
    #         else:
    #             raise ValueError(f"  Wrong options: {a}")
    #     else:
    #         unprocessed[key] = value
    #     width -= margin * 2
    #     height -= margin * 2

    # kwargs.clear()
    # kwargs |= unprocessed
    # logger.info(kwargs)
    # return Options(
    #     encode=encode,
    #     poly=poly,
    #     shadow=shadow,
    #     oxygen=oxygen,
    #     HB=HB,
    #     OH=OH,
    #     hydrogen=hydrogen,
    #     arrows=arrows,
    #     bgcolor=bgcolor,
    #     proj=proj,
    #     width=width,
    #     height=height,
    #     margin=margin,
    #     unprocessed=unprocessed,
    # )
    if options.shadow is True:
        options.shadow = "#8881"

    logger.info(f"options.rotate: {options.rotate}")
    for value in options.rotate:
        axis = value[0]
        angle = radians(float(value[1:]))
        cosx = cos(angle)
        sinx = sin(angle)
        if axis in "xX":
            R = np.array([[1, 0, 0], [0, cosx, sinx], [0, -sinx, cosx]])
        elif axis in "yY":
            R = np.array([[cosx, 0, -sinx], [0, 1, 0], [sinx, 0, cosx]])
        elif axis in "zZ":
            R = np.array([[cosx, sinx, 0], [-sinx, cosx, 0], [0, 0, 1]])
        else:
            raise ValueError(f"  Wrong options: {value}")
        options.proj = options.proj @ R

    if options.unprocessed:
        logger.warning(f"Unprocessed options: {options.unprocessed}")
        for key, value in options.unprocessed.items():
            logger.warning(f"  {key}: {value}")
    if options.width > 0:
        options.width = int(options.width)
    if options.height > 0:
        options.height = int(options.height)
    return options


def render_lattice_sites(genice: GenIce3, renderer: Render, options: dict):
    "A. Output molecular positions in PNG/SVG format."
    logger = getLogger()
    if options.hydrogen > 0 or options.arrows:
        # draw everything in hook6
        return
    offset = np.zeros(3)

    for i in range(3):
        options.proj[i] /= np.linalg.norm(options.proj[i])
    options.proj = np.linalg.inv(options.proj)

    cellmat = genice.cell
    projected = cellmat @ options.proj
    pos = genice.lattice_sites
    prims = []
    RO = options.oxygen  # nm
    RHB = options.oxygen * options.HB  # nm
    xmin, xmax, ymin, ymax = draw_cell(prims, projected)
    if options.poly:
        for ring in cycles_iter(genice.graph, 8, pos=pos):
            nedges = len(ring)
            deltas = np.zeros((nedges, 3))
            d2 = np.zeros(3)
            for k, i in enumerate(ring):
                d = pos[i] - pos[ring[0]]
                d -= np.floor(d + 0.5)
                deltas[k] = d
            comofs = np.sum(deltas, axis=0) / len(ring)
            deltas -= comofs
            com = pos[ring[0]] + comofs
            com -= np.floor(com)
            # rel to abs
            com = com @ projected
            deltas = deltas @ projected
            prims.append([com, "P", deltas, {"fillhs": hue_sat[nedges]}])  # line
    else:
        for i, j in genice.graph.edges():
            vi = pos[i]
            d = pos[j] - pos[i]
            d -= np.floor(d + 0.5)
            clipped = clip_cyl(vi @ projected, RO, (vi + d) @ projected, RO, RHB)
            if clipped is not None:
                prims.append(clipped + [RHB, {"fill": "#fff"}])  # line
            # If the bond intersects the cell boundary,
            if np.linalg.norm(vi + d - pos[j]) > 0.01:
                vj = pos[j]
                d = pos[i] - pos[j]
                d -= np.floor(d + 0.5)
                clipped = clip_cyl(vj @ projected, RO, (vj + d) @ projected, RO, RHB)
                if clipped is not None:
                    prims.append(clipped + [RHB, {"fill": "#fff"}])  # line
        for i, v in enumerate(pos):
            p = v @ projected
            prims.append([p, "C", RO, {}])  # circle
            # if the atom is on the boundary,
            for j in range(3):
                if p[j] == 0:
                    logger.info(f"On the boundary {j}")
                    q = v.copy()
                    q[j] = 1.0
                    p = q @ projected
                    prims.append([p, "C", RO, {}])  # circle
    # size of the object
    xsize = xmax - xmin
    ysize = ymax - ymin
    zoom = 100
    if options.width > 0:
        zoom = options.width / xsize
        if options.height > 0:
            z2 = options.height / ysize
            if z2 < zoom:
                zoom = z2
                xsize = options.width / zoom
                xcenter = (xmax + xmin) / 2
                xmin, xmax = xcenter - xsize / 2, xcenter + xsize / 2
            else:
                ysize = options.height / zoom
                ycenter = (ymax + ymin) / 2
                ymin, ymax = ycenter - ysize / 2, ycenter + ysize / 2
    elif options.height > 0:
        zoom = options.height / ysize
    logger.debug("Zoom {0} {1}x{2}".format(zoom, zoom * xsize, zoom * ysize))
    # margin in object scale
    rmargin = options.margin / zoom
    output = renderer(
        prims,
        RO,
        shadow=options.shadow,
        topleft=np.array((xmin - rmargin, ymin - rmargin)),
        size=(xsize + rmargin * 2, ysize + rmargin * 2),
        zoom=zoom,
        bgcolor=options.bgcolor,
        encode=options.encode,
    )
    return output


def render_atomic_sites(genice: GenIce3, renderer: Render, options: dict):
    "A. Output atomic positions in PNG/SVG format."
    logger = getLogger()

    filloxygen = {
        "stroke_width": 1,
        "stroke": "#444",
        "fill": "#f00",
        # "stroke_linejoin": "round",
        # "stroke_linecap" : "round",
        # "fill_opacity": 1.0,
    }
    fillhydrogen = {
        "stroke_width": 1,
        "stroke": "#444",
        "fill": "#0ff",
        # "stroke_linejoin": "round",
        # "stroke_linecap" : "round",
        # "fill_opacity": 1.0,
    }
    lineOH = {
        "stroke_width": 1,
        "stroke": "#444",
        "fill": "#fff",
    }
    lineHB = {
        "stroke_width": 1,
        "stroke": "#444",
        "fill": "#ff0",
    }
    arrow = {
        "stroke_width": 3,
        "stroke": "#fff",
    }
    offset = np.zeros(3)

    # Projection to the viewport
    for i in range(3):
        options.proj[i] /= np.linalg.norm(options.proj[i])
    options.proj = np.linalg.inv(options.proj)

    cellmat = genice.cell
    projected = cellmat @ options.proj

    prims = []
    RO = options.oxygen  # nm
    RHB = options.oxygen * options.HB  # nm
    ROH = options.oxygen * options.OH  # nm
    RH = options.oxygen * options.hydrogen  # nm
    xmin, xmax, ymin, ymax = draw_cell(prims, projected)
    if options.arrows:
        pos = genice.lattice_sites
        for i, j in genice.graph.edges():
            vi = pos[i]
            d = pos[j] - pos[i]
            d -= np.floor(d + 0.5)
            clipped = clip_cyl(
                vi @ projected, RO, (vi + d) @ projected, RO, 0.0
            )  # line
            if clipped is not None:
                prims.append(clipped + [0.0, {"stroke": "#fff"}])  # line
            if np.linalg.norm(vi + d - pos[j]) > 0.01:
                vj = pos[j]
                d = pos[i] - pos[j]
                d -= np.floor(d + 0.5)
                clipped = clip_cyl((vj + d) @ projected, RO, vj @ projected, RO, 0.0)
                if clipped is not None:
                    prims.append(clipped + [0.0, {"stroke": "#fff"}])  # line
        for i, v in enumerate(pos):
            prims.append([np.dot(v, projected), "C", RO, {}])  # circle
    else:
        water_model = parse_water_model_option("3site")
        waters = genice.water_molecules(water_model=water_model)

        # draw water molecules
        for water in waters.values():
            O = water.sites[0]
            H0 = water.sites[1]
            H1 = water.sites[2]
            prims.append([O @ options.proj, "C", RO, filloxygen])  # circle
            prims.append([H0 @ options.proj, "C", RH, fillhydrogen])  # circle
            prims.append([H1 @ options.proj, "C", RH, fillhydrogen])  # circle
            # clipped cylinder
            clipped = clip_cyl(O @ options.proj, RO, H0 @ options.proj, RH, ROH)
            if clipped is not None:
                prims.append(clipped + [ROH, lineOH])
            clipped = clip_cyl(O @ options.proj, RO, H1 @ options.proj, RH, ROH)
            if clipped is not None:
                prims.append(clipped + [ROH, lineOH])
        # draw HBs
        for i, j in genice.digraph.edges():
            if i in waters and j in waters:  # edge may connect to the dopant
                O = waters[j].sites[0]
                H0 = waters[i].sites[1]
                H1 = waters[i].sites[2]
                d0 = H0 - O
                d1 = H1 - O
                rr0 = d0 @ d0
                rr1 = d1 @ d1
                if rr0 < rr1 and rr0 < 0.245**2:
                    clipped = clip_cyl(O @ options.proj, RO, H0 @ options.proj, RH, RHB)
                    if clipped is not None:
                        prims.append(clipped + [RHB, lineHB])
                elif rr1 < rr0 and rr1 < 0.245**2:
                    clipped = clip_cyl(O @ options.proj, RO, H1 @ options.proj, RH, RHB)
                    if clipped is not None:
                        prims.append(clipped + [RHB, lineHB])
                # else:
                #     logger.debug(
                #         (np.linalg.norm(d["vector"]), rr0, rr1, 0.245**2)
                #     )
    xsize = xmax - xmin
    ysize = ymax - ymin
    zoom = 100
    rmargin = options.margin / zoom
    output = renderer(
        prims,
        RO,
        shadow=options.shadow,
        topleft=np.array((xmin - rmargin, ymin - rmargin)),
        size=(xsize + rmargin * 2, ysize + rmargin * 2),
        bgcolor=options.bgcolor,
        encode=options.encode,
    )
    return output


def dumps(genice: GenIce3, **kwargs):
    options = parse_options(**kwargs)
    renderer = Render
    if options.hydrogen > 0 or options.arrows:
        return render_atomic_sites(genice, renderer, options)
    else:
        return render_lattice_sites(genice, renderer, options)


def dump(genice: GenIce3, file: TextIOWrapper = sys.stdout, **kwargs):
    file.write(dumps(genice, **kwargs))

