#!/usr/bin/env bash
#
# setup_and_sync_repos.sh
#
# Run this ON YOUR MAC (not in the cloud sandbox), since it needs to
# read/write the local project folders and push to GitHub using your
# existing gh/git credentials.
#
# PHASE 1 - Two brand-new projects: git init, commit the existing files
#           as-is (nothing deleted/modified), push to the matching remote.
# PHASE 2 - Four already-initialized projects: commit any pending changes
#           and push.
#
# Requirements: git, and GitHub CLI `gh` already authenticated
#   (check: gh auth status ; if needed: gh auth login)
#
# Usage:
#   chmod +x setup_and_sync_repos.sh
#   ./setup_and_sync_repos.sh

set -uo pipefail

# ---------- Configuration ----------

NEW_REPOS=(
  "/Users/forbesa/Pathos_Projects/BD2.0/CD70_mAb|https://github.com/pathos-inc/p0138_BD2_CD70_ADC"
  "/Users/forbesa/Pathos_Projects/BD2.0/TROP2_HER3_ADC|https://github.com/pathos-inc/p0137_TROP2_HER3_ADC"
)

EXISTING_REPOS=(
  "/Users/forbesa/Pathos_Projects/p0130_MET_inhibitor"
  "/Users/forbesa/Pathos_Projects/p0122_poci_mm"
  "/Users/forbesa/Pathos_Projects/p0085_Rain_Diligence"
  "/Users/forbesa/Pathos_Projects/BD2.0/ER_PROTAC"
)

DEFAULT_BRANCH="main"
COMMIT_MSG_INIT="Initial commit of existing project files"

# ---------- Helpers ----------

log()  { printf '\n\033[1;34m==>\033[0m %s\n' "$1"; }
warn() { printf '\033[1;33m  ! %s\033[0m\n' "$1"; }
err()  { printf '\033[1;31m  x %s\033[0m\n' "$1"; }
ok()   { printf '\033[1;32m  + %s\033[0m\n' "$1"; }

check_requirements() {
  command -v git >/dev/null 2>&1 || { err "git not found. Install git first."; exit 1; }
  command -v gh  >/dev/null 2>&1 || { err "GitHub CLI (gh) not found. Install with: brew install gh"; exit 1; }
  gh auth status >/dev/null 2>&1 || { err "gh is not authenticated. Run: gh auth login"; exit 1; }
}

warn_large_files() {
  local dir="$1"
  local big
  big=$(find "$dir" -type f -size +90M -not -path '*/.git/*' 2>/dev/null)
  if [[ -n "$big" ]]; then
    warn "Files over 90MB found in $dir (GitHub hard-blocks anything over 100MB):"
    echo "$big" | sed 's/^/      /'
    warn "These may need Git LFS or a .gitignore exclusion, or the push below may fail."
  fi
}

init_and_push_new_repo() {
  local dir="$1" remote="$2"

  if [[ ! -d "$dir" ]]; then
    err "Directory not found: $dir - skipping."
    return 1
  fi

  log "Setting up new repo: $dir -> $remote"
  cd "$dir" || return 1

  warn_large_files "$dir"

  if [[ -d ".git" ]]; then
    warn "$dir already has a .git folder - skipping init, will just verify remote/commit/push."
  else
    git init -b "$DEFAULT_BRANCH"
    ok "git init done"
  fi

  local current_branch
  current_branch=$(git symbolic-ref --short -q HEAD || echo "")
  if [[ -z "$current_branch" ]]; then
    git checkout -b "$DEFAULT_BRANCH"
  fi

  if git remote get-url origin >/dev/null 2>&1; then
    local existing_remote
    existing_remote=$(git remote get-url origin)
    if [[ "$existing_remote" != "$remote" ]]; then
      warn "origin already points to $existing_remote - updating to $remote"
      git remote set-url origin "$remote"
    fi
  else
    git remote add origin "$remote"
    ok "remote 'origin' added"
  fi

  git add -A
  if git diff --cached --quiet; then
    ok "Nothing new to commit in $dir"
  else
    git commit -m "$COMMIT_MSG_INIT"
    ok "Committed existing files"
  fi

  if git ls-remote --exit-code --heads "$remote" "$DEFAULT_BRANCH" >/dev/null 2>&1; then
    log "Remote already has a '$DEFAULT_BRANCH' branch - fetching and merging before push"
    git fetch origin "$DEFAULT_BRANCH"
    if ! git merge --allow-unrelated-histories -m "Merge remote $DEFAULT_BRANCH into local" "origin/$DEFAULT_BRANCH"; then
      err "Merge conflict in $dir - resolve manually, then run:"
      err "  cd \"$dir\" && git push -u origin $DEFAULT_BRANCH"
      return 1
    fi
  fi

  git push -u origin "$DEFAULT_BRANCH"
  ok "Pushed $dir to $remote"
}

sync_existing_repo() {
  local dir="$1"

  if [[ ! -d "$dir" ]]; then
    err "Directory not found: $dir - skipping."
    return 1
  fi

  log "Checking for changes: $dir"
  cd "$dir" || return 1

  if [[ ! -d ".git" ]]; then
    err "$dir has no .git folder - expected an existing repo here. Skipping."
    return 1
  fi

  if ! git remote get-url origin >/dev/null 2>&1; then
    err "$dir has no 'origin' remote configured - skipping push."
    return 1
  fi

  warn_large_files "$dir"

  local branch
  branch=$(git symbolic-ref --short -q HEAD || echo "$DEFAULT_BRANCH")

  git pull --rebase origin "$branch" 2>/dev/null || warn "Could not pull/rebase (branch may not exist on remote yet, or there's a conflict) - continuing"

  if [[ -z "$(git status --porcelain)" ]]; then
    ok "No changes in $dir"
    return 0
  fi

  git add -A
  git commit -m "Sync changes - $(date +'%Y-%m-%d %H:%M %Z')"
  ok "Committed changes in $dir"

  git push origin "$branch"
  ok "Pushed $dir"
}

# ---------- Main ----------

check_requirements

log "PHASE 1: Initializing and pushing new repositories"
for entry in "${NEW_REPOS[@]}"; do
  IFS='|' read -r dir remote <<< "$entry"
  init_and_push_new_repo "$dir" "$remote" || err "Failed on $dir"
done

log "PHASE 2: Syncing existing repositories"
for dir in "${EXISTING_REPOS[@]}"; do
  sync_existing_repo "$dir" || err "Failed on $dir"
done

log "Done."
