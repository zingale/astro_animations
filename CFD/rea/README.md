This requires the `grid_plot.py` module from:

https://github.com/Open-Astrophysics-Bookshelf/numerical_exercises

To make the movies:

```
mkmovie.py -N 25  -o rea_nolimit `cat list_nolimit.txt | awk '{printf "%s ",$1}'`
```
