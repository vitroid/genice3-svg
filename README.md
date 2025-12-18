# genice3-svg

A GenIce3 exporter plugin for SVG and PNG visualization.

## Installation

```bash
pip install genice3-svg
```

## Usage

### SVG Export

```python
from genice3.plugin import UnitCell, Exporter
from genice3.genice import GenIce3

# Create a GenIce3 instance
genice = GenIce3(unitcell=UnitCell("ice1h_unit"))

# Export to SVG
with open("output.svg", "w") as f:
    Exporter("svg").dump(genice, file=f)
```

Or from command line:

```bash
genice3 CS2 -r 3 3 3 -e svg[rotatex=30,shadow] > CS2.svg
```

### PNG Export

```python
from genice3.plugin import UnitCell, Exporter
from genice3.genice import GenIce3

# Create a GenIce3 instance
genice = GenIce3(unitcell=UnitCell("ice1h_unit"))

# Export to PNG
with open("output.png", "wb") as f:
    Exporter("png").dump(genice, file=f)
```

Or from command line:

```bash
genice3 CS2 -r 3 3 3 -e png[shadow,bg=#f00] > CS2.png
```

## Options

- `rotatex=30`, `rotatey=30`, `rotatez=30`: Rotation angles
- `polygon`: Draw polygons instead of a ball and stick model
- `arrows`: Draw the hydrogen bonds with arrows
- `shadow=#8881`: Draw shadows behind balls
- `bg=#f00`: Specify the background color
- `O=0.06`: Size of the oxygen atom in nm
- `H=0`: Size of the hydrogen atom (relative to that of oxygen)
- `HB=0.4`: Radius of HB relative to that of oxygen
- `OH=0.5`: Radius of OH covalent bond relative to that of oxygen
- `width=0`, `height=0`: Image dimensions in pixels
- `margin=0`: Margin around the cell cube in pixels

## Requirements

- Python >= 3.9
- genice3 >= 3.0.a0
- svgwrite >= 1.4
- Pillow >= 10.0
- numpy >= 1.24
- cycless >= 0.4.2
- networkx >= 2.0

## License

MIT

