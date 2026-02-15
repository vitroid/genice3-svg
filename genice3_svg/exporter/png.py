# coding: utf-8
"""
GenIce format plugin to generate a PNG file.

Usage:
    % genice2 CS2 -r 3 3 3 -f png[shadow:bg=#f00] > CS2.png

Options:
    rotatex=30
    rotatey=30
    rotatez=30
    shadow         Draw shadows behind balls.
    bg=#f00        Specify the background color.
    H=0            Size of the hydrogen atom (relative to that of oxygen)
    O=0.06         Size of the oxygen atom in nm.
    HB=0.4         Radius of HB relative to that of oxygen
    OH=0.5         Radius of OH colvalent bond relative to that of oxygen
    width=0        (Pixel)
    height=0       (Pixel)
"""


from genice3_svg.exporter.svg import (
    parse_options,
    render_atomic_sites,
    render_lattice_sites,
)
from genice3_svg.render_png import Render
from genice3.genice import GenIce3
from io import TextIOWrapper
import sys

desc = {
    "ref": {},
    "brief": "PNG (Portable Network Graphics).",
    "usage": __doc__,
}




def dumps(genice: GenIce3, **kwargs):
    kwargs["encode"] = False
    options = parse_options(**kwargs)
    renderer = Render
    if options.hydrogen > 0 or options.arrows:
        return render_atomic_sites(genice, renderer, options)
    else:
        return render_lattice_sites(genice, renderer, options)


def dump(genice: GenIce3, file=sys.stdout, **kwargs):
    kwargs["encode"] = True
    # dumpsの返すImageをbinaryデータとして書き込みたい。
    image = dumps(genice, **kwargs)
    image.save(file, format="PNG")






