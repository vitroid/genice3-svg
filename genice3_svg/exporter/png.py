# coding: utf-8
"""
GenIce format plugin to generate a PNG file.

Usage:
    % genice3 CS2 -r 3 3 3 -e png[shadow,bgcolor=#f00] > CS2.png

Options:
    svg プラグインと同じオプション（rotate, polygon, arrows, shadow, bgcolor, O, H, HB, OH, width, height, margin）
"""

from genice3_svg.exporter.svg import (
    Options,
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
    """
    PNG 形式で出力。CLI からは parse_options の処理済み辞書が **kwargs で渡される。
    """
    options = Options(**kwargs, encode=False)
    renderer = Render
    if options.H > 0 or options.arrows:
        return render_atomic_sites(genice, renderer, options)
    else:
        return render_lattice_sites(genice, renderer, options)


def dump(genice: GenIce3, file: TextIOWrapper = sys.stdout, **kwargs):
    image = dumps(genice, **kwargs)
    image.save(file, format="PNG")






