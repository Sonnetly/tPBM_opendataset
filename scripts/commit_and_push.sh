#!/usr/bin/env bash
# Commit every non-ignored project change and push to the selected branch.
set -euo pipefail

if (( $# > 3 )); then
    printf 'Usage: bash scripts/commit_and_push.sh [message] [remote] [branch]\n' >&2
    exit 2
fi

message="${1:-Update project analysis files}"
remote="${2:-origin}"
branch="${3:-main}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(git -C "$script_dir" rev-parse --show-toplevel)"
cd -- "$repo_root"

git remote get-url "$remote" >/dev/null
git check-ref-format "refs/heads/$branch"

# Check the remote before staging or committing. Never discard local changes.
git fetch --no-tags "$remote" "refs/heads/$branch"
if ! git merge-base --is-ancestor FETCH_HEAD HEAD; then
    printf 'Remote %s/%s has commits missing locally. Integrate them, then retry. No files were staged or committed.\n' "$remote" "$branch" >&2
    exit 1
fi

git add -A -- .
git diff --cached --check
if git diff --cached --quiet; then
    printf 'No new file changes to commit. Pushing existing local commits.\n'
else
    git diff --cached --stat
    git commit -m "$message"
fi

# A concurrent remote change is rejected by normal fast-forward push rules.
git push "$remote" "HEAD:refs/heads/$branch"
local_commit="$(git rev-parse HEAD)"
remote_commit="$(git ls-remote --exit-code "$remote" "refs/heads/$branch" | awk '{print $1}')"
if [[ "$local_commit" != "$remote_commit" ]]; then
    printf 'Remote branch changed before verification. Expected %s; observed %s.\n' "$local_commit" "$remote_commit" >&2
    exit 1
fi
printf 'Verified %s/%s at %s\n' "$remote" "$branch" "$local_commit"
