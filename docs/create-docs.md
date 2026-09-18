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

```bash
cmake --build build-docs --target doc_lint
```

`doc_lint` checks every doc comment in the tree and publishes nothing. The two
targets exist because Doxygen validates a comment only while emitting that
entity's output, so one configuration cannot both publish a curated set and
check the rest: `Doxyfile`'s `FILE_PATTERNS` is what a user reads, and
`Doxyfile-lint` widens the input to the whole tree, turns on `EXTRACT_PRIVATE`,
and writes XML to a throwaway directory. It needs nothing compiled and takes
under a second. CI runs it on every pull request.

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

`.github/actions/build-site` runs exactly these steps with m.css pinned to a
commit. Every pull request builds the site and uploads it as the `website`
artifact; `docs.yml` builds it again on a `v*` tag and commits the result to the
`gh-pages` branch, which is what the published site serves.

@note m.css is stricter than Doxygen and aborts where Doxygen only warns — an
undocumented `friend` declaration, for instance, reaches its `parse_func` as a
member with an empty `argsstring` and raises
`AttributeError: 'NoneType' object has no attribute 'endswith'`. Wrap such
declarations in `@cond INTERNAL` / `@endcond`. The pull-request build is what
catches this before a tag does.
