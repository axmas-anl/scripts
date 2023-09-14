## CIF correlation plotting
A simple script for plotting some atomic correlation vectors in nexpy, likely for the purposes of superimposing on a PDF dataset.

## Usage
Define an AtomAnnotations object and us the `plot` method to add annotations to a plot, and `clear` to clear them. 
```
s = read_cif("/path/to/a/crystallographic/information/file.cif")
a = AtomAnnotations(s, plotview.ax)
a.plot('xy0', 0.0)
```

### The plot method
The plot method takes two arguments, a strin reprenting the plane and an offset. 
The string should be one of "xy0", "x0z", or "0xz". 
The offset should be the _actual_ value in the vertical direction, i.e. to plot correlations close to the [x, 0.575, z] plane then use `plot("x0z", 0.575)`.
Optionally one can provide a threshold, i.e. how close to the plane the correlation needs to be to appear. By default plot will add any correlation within 0.5 Angstroms. 

