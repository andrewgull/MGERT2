# MGERT2

[![Unit Tests](https://github.com/andrewgull/MGERT2/actions/workflows/units.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/units.yml)
[![Code Coverage](https://github.com/andrewgull/MGERT2/actions/workflows/test-and-coverage.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/test-and-coverage.yml)
[![Small integration test](https://github.com/andrewgull/MGERT2/actions/workflows/integration.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/integration.yml)
[![Snakemake Dry Run Check](https://github.com/andrewgull/MGERT2/actions/workflows/snakemake-dry-run.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/snakemake-dry-run.yml)
[![Snakefile Formatting Check](https://github.com/andrewgull/MGERT2/actions/workflows/snakefmt.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/snakefmt.yml)
[![Python Code Style](https://github.com/andrewgull/MGERT2/actions/workflows/pycodestyle.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/pycodestyle.yml)
[![Markdown Style](https://github.com/andrewgull/MGERT2/actions/workflows/markdown-lint.yml/badge.svg)](https://github.com/andrewgull/MGERT2/actions/workflows/markdown-lint.yml)
![CodeRabbit Pull Request Reviews](https://img.shields.io/coderabbit/prs/github/andrewgull/MGERT2?utm_source=oss&utm_medium=github&utm_campaign=andrewgull%2FMGERT2&labelColor=171717&color=FF570A&link=https%3A%2F%2Fcoderabbit.ai&label=CodeRabbit+Reviews)

## Introduction

This is a rework of [MGERT](https://github.com/andrewgull/MGERT) -
aka a big-ass 6 year old python script which was slain by "dependency hell".

Check the Issues tab and [Project requirements](docs/requirements.md) to get started.

Remember that `dev` and `main` branches are protected - do not push changes to them.
Create a new branch to contribute to the project.

Follow [this guide](https://gist.github.com/joshbuchea/6f47e86d2510bce28f8e7f42ae84c716)
to write better commit messages.

## What you need to install

To ensure we all use exactly the same dependencies to run the workflow,
tests etc., we use [Pixi](https://pixi.sh/).

Intall it on your laptop using this [guide](https://pixi.prefix.dev/latest/installation/).

## How to run

Having done this, you'll be able to run all the commands listed below:

```bash
# format Snakefile
# check for the visual style
pixi run format-snakefile

# format Python scripts
pixi run format-scripts

# lint Snakefile:
# check for logic and best practices
pixi run lint

# rulegraph (shows general workflow structure)
# to images/
pixi run rulegraph

# run unit tests
pixi run tests

# snakemake dry run
pixi run dry-run

# snakemake full run (4 cores)
pixi run full-run
```

To get full list of available tasks, run:

```bash
pixi task list
```

Also, see [pixi.toml](pixi.toml) for more details.
