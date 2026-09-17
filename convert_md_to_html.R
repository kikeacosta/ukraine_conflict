# ==============================================================================
# Convert Markdown Reports to HTML (using commonmark - a real CommonMark
# parser - instead of hand-rolled regex, which was never actually run and
# would have mishandled nested lists, tables, and code blocks)
# ==============================================================================

suppressMessages(library(commonmark))

md_files <- c(
  "IMPLEMENTATION_STATUS.md",
  "MULTISTATE_ANALYSIS_README.md",
  "MISSING_COMBATANT_IMPUTATION_ANALYSIS.md",
  "ualosses_multistate_methodology.md",
  "results_guide.md",
  "MULTISTATE_WORKFLOW.md",
  "PHASE_1_VALIDATION_GUIDE.md",
  "project_summary.md",
  "MULTISTATE_IMPROVEMENTS_ROADMAP.md",
  "QUICK_IMPROVEMENTS_SUMMARY.md",
  "HTML_PORTAL_README.md",
  "README.md",
  "SESSION_HANDOFF_2026-09-17.md"
)

template <- '<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>%s - Ukraine Conflict Analysis</title>
<style>
:root {
  --primary: #2c3e50; --secondary: #3498db; --accent: #e74c3c;
  --light: #ecf0f1; --dark: #34495e; --text: #2c3e50; --bg: #ffffff;
}
@media (prefers-color-scheme: dark) {
  :root { --text: #ecf0f1; --bg: #1a1a1a; --light: #2a2a2a; }
}
* { box-sizing: border-box; }
body {
  font-family: "Segoe UI", Tahoma, Geneva, Verdana, sans-serif;
  line-height: 1.7; color: var(--text); background-color: var(--bg);
  padding: 20px; margin: 0;
}
.container { max-width: 900px; margin: 0 auto; }
header {
  background: linear-gradient(135deg, var(--primary) 0%%, var(--secondary) 100%%);
  color: white; padding: 30px; border-radius: 8px; margin-bottom: 30px; text-align: center;
}
header h1 { font-size: 1.8em; margin: 0 0 8px 0; }
header p { opacity: 0.9; margin: 0; }
.nav-bar {
  margin-bottom: 30px; padding: 12px; background-color: var(--light);
  border-radius: 5px; text-align: center;
}
.nav-bar a { color: var(--secondary); text-decoration: none; margin: 0 12px; font-weight: bold; }
.nav-bar a:hover { text-decoration: underline; }
main h1 { color: var(--primary); margin: 30px 0 15px; font-size: 1.6em; padding-bottom: 8px; border-bottom: 2px solid var(--secondary); }
main h2 { color: var(--secondary); margin: 25px 0 12px; font-size: 1.35em; }
main h3 { color: var(--dark); margin: 20px 0 10px; font-size: 1.15em; }
main h4 { color: var(--dark); margin: 16px 0 8px; font-size: 1.0em; }
p { margin: 0 0 14px; }
ul, ol { margin: 0 0 14px 28px; padding: 0; }
li { margin-bottom: 6px; }
code { background-color: var(--light); padding: 2px 6px; border-radius: 3px; font-family: "Consolas","Courier New",monospace; font-size: 0.9em; }
pre { background-color: var(--light); padding: 14px; border-radius: 5px; overflow-x: auto; margin: 14px 0; border-left: 4px solid var(--secondary); }
pre code { background: none; padding: 0; }
table { width: 100%%; border-collapse: collapse; margin: 18px 0; display: block; overflow-x: auto; }
th, td { padding: 10px 12px; text-align: left; border-bottom: 1px solid var(--light); }
th { background-color: var(--secondary); color: white; font-weight: bold; }
tr:hover { background-color: var(--light); }
blockquote { border-left: 4px solid var(--secondary); padding: 4px 15px; margin: 14px 0; color: var(--dark); background: var(--light); }
strong { font-weight: 600; }
a { color: var(--secondary); }
hr { border: none; height: 2px; background-color: var(--light); margin: 28px 0; }
footer { text-align: center; margin-top: 50px; padding-top: 20px; border-top: 2px solid var(--light); color: var(--dark); font-size: 0.85em; }
</style>
</head>
<body>
<div class="container">
<header><h1>%s</h1><p>Ukraine Conflict Mortality Analysis - Documentation</p></header>
<div class="nav-bar"><a href="index.html">&larr; Back to Index</a></div>
<main>
%s
</main>
<footer><p>Ukraine Conflict Mortality Analysis | Documentation | <a href="index.html">All Reports</a></p></footer>
</div>
</body>
</html>'

convert_one <- function(md_file) {
  src <- file.path("documents", md_file)
  if (!file.exists(src)) {
    message("MISSING: ", md_file)
    return(invisible(NULL))
  }
  md_text <- paste(readLines(src, warn = FALSE), collapse = "\n")

  # Strip the title's leading "# " line from the body since it's rendered
  # separately in the <header>, to avoid a duplicate title at the top.
  first_line <- sub("^#\\s+", "", strsplit(md_text, "\n")[[1]][1])
  title <- if (grepl("^#\\s", md_text)) first_line else tools::file_path_sans_ext(md_file)

  body_html <- markdown_html(md_text, extensions = c("table", "strikethrough", "autolink"))

  out_name <- paste0(tools::file_path_sans_ext(md_file), ".html")
  out_path <- file.path("documents", out_name)
  writeLines(sprintf(template, title, title, body_html), out_path, useBytes = TRUE)
  message(sprintf("OK: %s -> %s (%d KB)", md_file, out_name, round(file.info(out_path)$size / 1024)))
}

for (f in md_files) convert_one(f)

message("\nDone.")
