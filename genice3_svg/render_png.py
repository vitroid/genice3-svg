from logging import getLogger
import sys
import io

import numpy as np
import PIL.ImageDraw as ImageDraw
import PIL.Image as Image


def cylinder(draw, v1_, v2_, r, **options):
    """
    draw a 3D cylinder
    """
    options = {"fill": "#fff", **options}
    draw.line(
        [int(x) for x in [v1_[0], v1_[1], v2_[0], v2_[1]]],
        width=int(r * 2),
        fill=options["fill"],
    )


def Render(
    prims,
    Rsphere,
    shadow=None,
    topleft=np.array([-1.0, -1.0]),
    size=(50, 50),
    zoom=200,
    vertices=None,
    vecs=None,
    bgcolor="#fff",
    encode=True,
):
    """
    Renter the image in PNG format and output to the stdout.
    Returns nothing.

    When vertex list is given, the coords in prims are not indicated in vectors but in indices
    Vecs are vectors not needed to be sorted (used in "L" command)
    """
    logger = getLogger()
    size = tuple((int(x * zoom + 0.5) for x in size))
    image = Image.new("RGB", size, bgcolor)
    draw = ImageDraw.Draw(image, "RGBA")
    # special treatment for post-K project
    # draw.rectangle([0,0,size[0]/2,size[1]], fill="#EF5FA7")
    # draw.rectangle([size[0]/2,0,size[0],size[1]], fill="#00A2FF")

    TL0 = np.zeros(3)
    TL0[:2] = topleft
    linedefaults = {
        "stroke_width": 2,
        "stroke": "#000",
        # "stroke_linejoin": "round",
        # "stroke_linecap" : "round",
    }
    filldefaults = {
        "stroke_width": 1,
        "stroke": "#000",
        "fill": "#0ff",
        # "stroke_linejoin": "round",
        # "stroke_linecap" : "round",
        # "fill_opacity": 1.0,
    }
    shadowdefaults = {
        "stroke_width": 0,
        # "fill": "#8881",
        "fill": shadow,
        # "fill_opacity": 0.08,
    }
    if shadow is not None:
        r = Rsphere
        Z = np.array([0, 0, 1.0])
        prims += [
            [prim[0] - Z * r * 1.4**j, prim[1] + "S", r * 1.4**j] + prim[3:]
            for j in range(1, 5)
            for prim in prims
            if prim[1] == "C"
        ]
    prims.sort(key=lambda x: -x[0][2])
    while len(prims) > 0:
        prim = prims.pop()
        if not (
            (-0.5 < prim[0][0] - topleft[0] < size[0] + 0.5)
            and (-0.5 < prim[0][1] - topleft[1] < size[1] + 0.5)
        ):
            continue
        if prim[1] == "L":
            if prim[4] == 0:
                options = {**linedefaults, **prim[5]}
                s = (prim[2][:2] - topleft) * zoom
                e = (prim[3][:2] - topleft) * zoom
                draw.line(
                    [int(s[0]), int(s[1]), int(e[0]), int(e[1])],
                    fill=options["stroke"],
                    width=options["stroke_width"],
                )
            else:
                options = {**filldefaults, **prim[5]}
                cylinder(
                    draw,
                    (prim[2] - TL0) * zoom,
                    (prim[3] - TL0) * zoom,
                    prim[4] * zoom,
                    **options
                )
        elif prim[1] == "L2":
            # new, simpler expression.
            # half relative vector is given
            if prim[3] == 0:
                options = {**linedefaults, **prim[4]}
                s = ((prim[0] + prim[2])[:2] - topleft) * zoom
                e = ((prim[0] - prim[2])[:2] - topleft) * zoom
                draw.line(
                    [int(s[0]), int(s[1]), int(e[0]), int(e[1])],
                    fill=options["stroke"],
                    width=options["stroke_width"],
                )
            else:
                options = {**filldefaults, **prim[4]}
                cylinder(
                    draw,
                    (prim[0] + prim[2] - TL0) * zoom,
                    (prim[0] - prim[2] - TL0) * zoom,
                    prim[3] * zoom,
                    **options
                )
        elif prim[1] == "C":
            options = {**filldefaults, **prim[3]}
            Rsphere = prim[2]
            center = (prim[0][:2] - topleft) * zoom
            r = Rsphere * zoom
            tl = center - r
            br = center + r
            draw.ellipse(
                [int(x) for x in [tl[0], tl[1], br[0], br[1]]], fill=options["fill"]
            )
        elif prim[1] == "CS":
            Rsphere = prim[2]
            options = {
                **shadowdefaults,
            }  # **prim[3] }
            # logger.info("{0}".format(options))
            center = (prim[0][:2] - topleft) * zoom
            r = Rsphere * zoom
            tl = center - r
            br = center + r

            draw.ellipse(
                [int(x) for x in [tl[0], tl[1], br[0], br[1]]], fill=options["fill"]
            )
    if encode:
        imgByteArr = io.BytesIO()
        image.save(imgByteArr, format="PNG")
        imgByteArr = imgByteArr.getvalue()
        return imgByteArr
    else:
        return image

