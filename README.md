# Research 

You can download this repository:

```bash
git clone https://github.com/scientific-analytics/researchpaper_climatevar.git
```

```bash
python -m venv ven
source venv/bin/activate
pip install -r requirements.txt
```

## Reproducible Research

Code folder: a Jupyter Notebook with the code to reproduce the figures

Manuscript: LateX code, pdf and word documents

## Tools

[Go here for VS Code and LaTeX installation](https://mathjiajia.github.io/vscode-and-latex/).

To convert your LaTeX document into Word document, [install Pandoc](https://pandoc.org/installing.html). Then, in the folder where your `.tex` file is located, run:

```bash
pandoc -s input.tex -o output.docx
```

You may need to ask Scientific Beta Helpdesk for authorizing `pandoc.exe` to run the first time.


# Contributing

## Issue

When creating an issue, please follow these guidelines for the title:

- **Bug report**: If you are reporting a bug, start the title with `[BUG]`. If the issue concerns the Data or the Front-End and not directly the library itself, add `DATA` or `FRONT-END` after the enclosed `[BUG]`.

- **Feature request**: If you are requesting a new feature, start the title with `[FEATURE]`. If the issue concerns the Data or the Front-End and not directly the library itself, add `DATA` or `FRONT-END` after the enclosed `[FEATURE]`.

The issue creation will add a ticket into the `ESG` Jira backlog.

## Process

### Bugfixes

- create a branch from the issue page
- branch naming convention: `bugfix/<issue number>` or `bugfix/<bugfix-name>`
- pull request to: `main`

### Enhancements

- create a branch from the issue page
- branch naming convention: `feature/<issue number>` or `feature/<enhancement-name>`
- pull request to: `main`