# Research Project Agent Instructions

## Project Philosophy

This repository is intended to support reproducible computational research.

The primary goals are:

1. Scientific accuracy
2. Reproducibility
3. Maintainability
4. Education
5. Extensibility

The user wants to become a better programmer and computational scientist.
Teaching is equally as important as producing working code.

---

# Collaboration Style

Always:

- Explain the plan before implementing it.
- Explain scientific reasoning.
- Explain programming concepts.
- Explain tradeoffs when multiple solutions exist.

Never:

- Make scientific assumptions without explicitly stating them.
- Rewrite large portions of the repository unless requested.
- Remove documentation.
- Modify raw data.

---

# Repository Organization

## notebooks/

Purpose:
Scientific workflow and narrative.

Keep notebooks readable.

Notebook cells should primarily:

- load data
- call functions
- display figures
- explain results

Avoid placing large amounts of implementation code inside notebooks.

---

## src/

Contains reusable source code.

Prefer:

- small functions
- descriptive names
- docstrings
- type hints where appropriate

---

## data/raw/

Never modify.

Treat as immutable.

---

## data/processed/

Generated from raw data.

May be regenerated at any time.

---

## results/

Generated outputs.

Do not edit manually.

---

## docs/

Keep documentation updated whenever significant changes occur.

---

# Coding Standards

Prefer:

- readability
- simplicity
- maintainability

Avoid:

- unnecessary optimization
- duplicated code
- hidden constants
- magic numbers

Document algorithms.

---

# Git Workflow

Prefer small commits.

Recommend commit messages.

Do not rewrite Git history.

Do not force push unless explicitly requested.

---

# Testing

Whenever practical:

- add tests
- validate outputs
- compare against expected values
- identify edge cases

---

# Documentation

Whenever functionality changes:

Update:

- README.md
- replication_journal.md
- assumptions.md

if applicable.

---

# Scientific Reproducibility

If a published paper is ambiguous:

1. Identify the ambiguity.
2. Explain possible interpretations.
3. Recommend one.
4. Document the decision.

Never silently change published methods.

---

# Learning

The user is learning.

Whenever possible:

- explain concepts
- explain code
- suggest improvements
- identify opportunities for refactoring
- recommend best practices

Do not simply generate code without explanation unless explicitly requested.

---

# Preferred Workflow

Understand

↓

Plan

↓

Discuss

↓

Implement

↓

Test

↓

Document

↓

Commit