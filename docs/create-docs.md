# Building the documentation

The pages under `docs/` are Doxygen input; every code example in them is a
`\snippet` of a compiled source under `docs/snippets/`, so the examples cannot
drift from the API.

## Locally

```bash
cmake -S . -B build-docs -DVINECOPULIB_BUILD_DOC=ON -DBUILD_TESTING=OFF
cmake --build build-docs --target doc
```

`doc` depends on `doc_snippets`, so the examples are compiled first, and Doxygen
runs with `WARN_AS_ERROR = YES` — an undocumented parameter or a broken
reference fails the build. CI runs the same target, and also runs each snippet
binary, on every pull request.

## The website

The published site uses the [m.css](https://github.com/mosra/m.css) theme rather
than Doxygen's own HTML. `Doxyfile-mcss.in` is configured into the build tree
next to `Doxyfile`, so:

```bash
python3 -m pip install jinja2 Pygments
sudo apt-get install texlive-latex-recommended texlive-latex-extra \
                     texlive-fonts-recommended texlive-fonts-extra \
                     preview-latex-style dvisvgm
git clone https://github.com/mosra/m.css
cd build-docs && python3 ../m.css/documentation/doxygen.py Doxyfile-mcss
```

The result lands in `build-docs/docs/` (`HTML_OUTPUT` is `.`). LaTeX is needed
because m.css renders `\f$...\f$` math to SVG.

@note m.css requires a Doxygen version it agrees with: with Doxygen 1.9.8 its
`parse_func` aborts on a `memberdef` that carries no `argsstring`
(`AttributeError: 'NoneType' object has no attribute 'endswith'`). Deployment is
therefore still manual; automating it needs a pinned Doxygen or a patched m.css.
