# libparanumal — Claude Code project instructions

## Stacked PRs

Always use `gh stack` when creating stacked PRs for this project.

Workflow:
1. From `main`: create the first feature branch and initialize the stack
   ```
   git checkout -b feat/<slug>
   gh stack init --adopt feat/<slug>
   ```
2. For each subsequent PR in the stack:
   ```
   gh stack add feat/<next-slug>
   # implement changes, then commit
   ```
3. Push and submit all PRs at once:
   ```
   gh stack push
   gh stack submit
   ```
4. When rebasing after upstream changes:
   ```
   gh stack rebase
   ```

PR branches follow the naming convention: `feat/acoustics-<slug>`
