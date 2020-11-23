# 1. Clone your fork:
git clone git@github.com:gauravj49/nf-core/chipseq.git
# 2. Add remote from original repository in your forked repository:
cd /home/rad/users/gaurav/projects/workflows/nfchipseq
git remote add upstream https://github.com/nf-core/chipseq.git
git fetch upstream
# 3. Updating your fork from original repo to keep up with their changes:
git pull upstream master


# Fork repository from the github
# https://github.com/nf-core/chipseq

# Set config values
git config user.name  "Gaurav Jain"
git config user.email "gauravj49@gmail.com"

# Check the config list
git config --list

# Get the status of the project and repository
git status

# Ignore files that should not go into the repository
# emacs -nw .gitignore
cat > .gitignore
annotation
docs
input
output
.gitignore
scripts/usage/00_git_setup_usage.sh

# Once you have done that, git now knows about your remote repository. You can then tell it to push (which is "upload") your commited files:
git push


