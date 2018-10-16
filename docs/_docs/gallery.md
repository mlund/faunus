---
layout: splash
toc: false
permalink: /gallery/

example_row1:
  - image_path: https://www.dropbox.com/s/spsfqnphqgutwrd/water.png?raw=1
    url: "https://github.com/mlund/faunus/blob/master/examples/water.yml"
    alt: "SPC/E water in the NPT ensemble"
    title: "SPC/E water"
    btn_label: "Input file"
    btn_class: "btn--inverse"
    excerpt: "Simple MC example of explicit water with pressure coupling"
  - image_path: https://www.dropbox.com/s/7ymnpf9t4w1nt48/membrane-3bead.jpg?raw=1
    url: "{{site.github.repository_url}}/blob/master/examples/membrane.yml"
    title: "Coarse Grained Bilayer"
    excerpt: "Three-bead lipid bilayer model spanning a periodic simulation box to form a bilayer. Wang-Landau sampling of bending modulus." 
---
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Gallery

This is a display of systems Faunus can be used to simulate. Hover the mouse over an
image to see more information and click to see the input file.

{% include example_row1 id="example_row1" %}

