(TeX-add-style-hook
 "simulations"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "standalone"
    "standalone10"
    "tikz"))
 :latex)

