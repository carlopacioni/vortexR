# Implementation
<!--
For non-trivial PRs, explain your implementation decisions to the reviewers.
Aim to write commit messages such that they provide all low-level detail, include them here if useful.

Every PR should reference a GitHub issue which contains a thorough description of the feature or bug.

Suggested format:

#<issue>
* Explanation
* Details pasted from commit messages where appropriate.
* If the issus is sufficiently addressed, add the "Close #<issue>" statement to auto-close issue on PR merge.
-->

## Pre-submission checklist
<!--
For non-trivial PRS, we encourage to send draft PRs for early feedback.
Before you submit the PR, follow the steps below to make review easier and faster.

If the PR has many small commits with uninformative messages, feel free
to squash them into one commit with a well written commit message
in a new branch and start a fresh PR:

git checkout -b <issue-number>-short-title
git reset --soft $(git merge-base master HEAD)
git add -A
git commit -m "Great commit message"
git push

-->
- [ ] The submitted branch has a meaningful name: `<GH issue ID>-<short-title>`
- [ ] **Commit messages** follow [best practices](https://cbea.ms/git-commit/).
- [ ] Add valid Roxygen docstrings to functions.
- [ ] Add **tests** to cover expected behaviour, all possible input argument options, edge cases.
      Add at least documented test stubs to cover new functionality and create a new issue as reminder.
- [ ] Have **pre-commit** installed locally or run pre-commit checks prior to committing with
      `pre-commit run --all-files`
- [ ] All checks on the PR have succeeded.
