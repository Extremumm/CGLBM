# -----------------------------
# Options affecting formatting.
# -----------------------------
with section("format"):
    # How wide to allow formatted cmake files
    line_width = 100

    # How many spaces to tab for indent
    tab_size = 4

    # If a positional argument group contains more than this many arguments, then
    # force it to a vertical layout.
    max_pargs_hwrap = 4

    # If a statement is wrapped to more than one line, than dangle the closing
    # parenthesis on its own line.
    dangle_parens = True

    # Format command names consistently as 'lower' or 'upper' case
    command_case = "lower"

    # Format keywords consistently as 'lower' or 'upper' case
    keyword_case = "upper"

    # If true, the parsers may infer whether or not an argument list is sortable
    # (without annotation).
    autosort = True

# ----------------------------
# Options affecting the linter
# ----------------------------
with section("lint"):
    # a list of lint codes to disable
    disabled_codes = ["C0103"]
