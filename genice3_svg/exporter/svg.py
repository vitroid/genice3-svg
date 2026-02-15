# coding: utf-8
"""
GenIce exporter plugin to generate a SVG file.

Usage:
    % genice3 CS2 -r 3 3 3 -e svg[rotatex=30,shadow] > CS2.svg

Options:
    --rotate y5 x5       rotation axes+angles (CLI: space or comma). API: list ["y5","x5"]
    --polygon        Draw polygons instead of a ball and stick model.
    --arrows         Draw the hydrogen bonds with arrows.
    --shadow #8881   Draw shadows behind balls.
    --bgcolor #f00   Specify the background color.
    --O 0.06         Radius of the oxygen atom
    --H 0            Size of the hydrogen atom (relative to that of oxygen)
    --HB 0.4         Radius of HB relative to that of oxygem
    --OH 0.5         Radius of OH colvalent bond relative to that of oxygem
    --width 0        (Pixel)
    --height 0       (Pixel)
    --margin 0       (Pixel)
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
from dataclasses import dataclass, field
from genice3.genice import GenIce3
import sys
from io import TextIOWrapper
from genice3.exporter import (
    parse_water_model_option,
)
from genice3.cli.pool_parser import (
    OptionDef,
    parse_options_generic,
    OPTION_TYPE_STRING,
    OPTION_TYPE_TUPLE,
    OPTION_TYPE_FLAG,
)
from typing import Dict, Any, Tuple

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


def _normalize_rotate(value) -> list:
    """
    rotate オプションを "y5", "x5" のリストに正規化する。
    CLI: --rotate y5 x5 または --rotate y5,x5。API: リスト ["y5","x5"]。
    """
    if isinstance(value, str):
        return [s.strip() for s in re.split(r"[\s,]+", value) if s.strip()]
    if isinstance(value, (list, tuple)):
        out = []
        for item in value:
            s = item if isinstance(item, str) else str(item)
            for part in re.split(r"[\s,]+", s):
                if part.strip():
                    out.append(part.strip())
        return out
    return [str(value)]


def rotation_processor(x) -> np.ndarray:
    """axis+angle のリスト (e.g. ["y5", "x5"]) から回転行列を生成。"""
    items = _normalize_rotate(x)
    mat = np.eye(3)
    for value in items:
        if not value:
            continue
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
            assert False, "  Wrong options."
        mat = mat @ R
    return mat


# svg プラグインが受け取るオプション定義。追加・削除はここだけ行えばよい。
SVG_OPTION_DEFS = (
    OptionDef("rotate", parse_type=OPTION_TYPE_TUPLE),
    OptionDef("polygon", parse_type=OPTION_TYPE_FLAG),
    OptionDef("arrows", parse_type=OPTION_TYPE_FLAG),
    OptionDef("shadow", parse_type=OPTION_TYPE_STRING),
    OptionDef("bgcolor", parse_type=OPTION_TYPE_STRING),
    OptionDef("O", parse_type=OPTION_TYPE_STRING),
    OptionDef("H", parse_type=OPTION_TYPE_STRING),
    OptionDef("HB", parse_type=OPTION_TYPE_STRING),
    OptionDef("OH", parse_type=OPTION_TYPE_STRING),
    OptionDef("width", parse_type=OPTION_TYPE_STRING),
    OptionDef("height", parse_type=OPTION_TYPE_STRING),
    OptionDef("margin", parse_type=OPTION_TYPE_STRING),
)


def parse_options(options: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    svgプラグインのオプションを処理

    --rotate y5 x5       rotation (CLI: y5 x5 or y5,x5). API: list ["y5","x5"]
    --polygon        Draw polygons instead of a ball and stick model.
    --arrows         Draw the hydrogen bonds with arrows.
    --shadow #8881   Draw shadows behind balls.
    --bgcolor #f00   Specify the background color.
    --O 0.06         Radius of the oxygen atom
    --H 0            Size of the hydrogen atom (relative to that of oxygen)
    --HB 0.4         Radius of HB (relative to that of oxygem)
    --OH 0.5         Radius of OH colvalent bond (relative to that of oxygem)
    --width 0        (Pixel)
    --height 0       (Pixel)
    --margin 0       (Pixel)

    Returns:
        (処理したオプション, 処理しなかったオプション)。未処理は次のプラグインへ。
    """
    options = dict(options)
    option_specs = {
        d.name: d.parse_type
        for d in SVG_OPTION_DEFS
        if d.parse_type is not None
    }
    post_processors = {
        "shadow": lambda x: "#8881" if x is True else x,
        "width": lambda x: int(x),
        "height": lambda x: int(x),
        "margin": lambda x: int(x),
        "O": lambda x: float(x),
        "H": lambda x: float(x),
        "HB": lambda x: float(x),
        "OH": lambda x: float(x),
        "rotate": lambda x: rotation_processor(x),
    }
    return parse_options_generic(options, option_specs, post_processors)



def _default_rotate():
    """parse_options で rotate が渡されない場合のデフォルト回転行列。"""
    return rotation_processor(("x0",))


@dataclass
class Options:
    rotate: np.ndarray = field(default_factory=_default_rotate)
    O: float = 0.06
    H: float = 0   
    HB: float = 0.4
    OH: float = 0.5
    width: int = 0
    height: int = 0
    margin: int = 0
    encode: bool = True
    command_line:str = ""
    polygon: bool = False
    arrows: bool = False
    shadow: str = "#8881"
    bgcolor: str = "#ff0000"

def render_lattice_sites(genice: GenIce3, renderer: Render, options: Options):
    "A. Output molecular positions in PNG/SVG format."
    logger = getLogger()
    if options.H > 0 or options.arrows:
        # draw everything in hook6
        return
    offset = np.zeros(3)

    for i in range(3):
        options.rotate[i] /= np.linalg.norm(options.rotate[i])
    options.rotate = np.linalg.inv(options.rotate)

    cellmat = genice.cell
    projected = cellmat @ options.rotate
    pos = genice.lattice_sites
    prims = []
    RO = options.O  # nm
    RHB = options.O * options.HB  # nm
    xmin, xmax, ymin, ymax = draw_cell(prims, projected)
    if options.polygon:
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
        options.rotate[i] /= np.linalg.norm(options.rotate[i])
    options.rotate = np.linalg.inv(options.rotate)

    cellmat = genice.cell
    projected = cellmat @ options.rotate

    prims = []
    RO = options.O  # nm
    RHB = options.O * options.HB  # nm
    ROH = options.O * options.OH  # nm
    RH = options.O * options.H  # nm
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
            prims.append([O @ options.rotate, "C", RO, filloxygen])  # circle
            prims.append([H0 @ options.rotate, "C", RH, fillhydrogen])  # circle
            prims.append([H1 @ options.rotate, "C", RH, fillhydrogen])  # circle
            # clipped cylinder
            clipped = clip_cyl(O @ options.rotate, RO, H0 @ options.rotate, RH, ROH)
            if clipped is not None:
                prims.append(clipped + [ROH, lineOH])
            clipped = clip_cyl(O @ options.rotate, RO, H1 @ options.rotate, RH, ROH)
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
                    clipped = clip_cyl(O @ options.rotate, RO, H0 @ options.rotate, RH, RHB)
                    if clipped is not None:
                        prims.append(clipped + [RHB, lineHB])
                elif rr1 < rr0 and rr1 < 0.245**2:
                    clipped = clip_cyl(O @ options.rotate, RO, H1 @ options.rotate, RH, RHB)
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
    options = Options(**kwargs)
    renderer = Render
    if options.H > 0 or options.arrows:
        return render_atomic_sites(genice, renderer, options)
    else:
        return render_lattice_sites(genice, renderer, options)


def dump(genice: GenIce3, file: TextIOWrapper = sys.stdout, **kwargs):
    file.write(dumps(genice, **kwargs))






