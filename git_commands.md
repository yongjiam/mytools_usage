## github push steps:
0. create a new repository in github account

## create ssh keys pairs and use it instead of http
1. check all ssh keys
ls -al ~/.ssh

2. generate ssh key pairs
ssh-keygen -t rsa -b 4096 -C "your_email@example.com"  ## path/to/ssh/keyname

3. add public key to github account
https://docs.github.com/en/github-ae@latest/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account
> click your profile photo, then click Settings \
> In the "Access" section of the sidebar, click  SSH and GPG keys \
> Click New SSH key or Add SSH key \
> In the "Title" field, add a descriptive label for the new key \
> Select the type of key, either authentication or signing \
> In the "Key" field, paste your public key \
> cat ~/.ssh/id_rsa.pub \
> Then select and copy the contents of the id_ed25519.pub file

5. add local key to ssh agent
ssh-add /path/to/your/private_key

start ssh agent if got error with ssh-add
eval "$(ssh-agent)"

5. check remote url: http or ssh
git remote -v

6. set remote origin curl as ssh
git remote set-url origin git@github.com:your_username/your_repository_name.git

7. git clone remote repository

8. cp or make changes to local git directory

9. initiate local git
git init

10. to stage the changes
git add .

11. registor changes
git commit -m "updates"

12. git push -u origin branch_name(main or master) ## -u flag in the push command (git push -u origin main) sets the upstream branch, but it's only necessary for the first push
or just: git push master/main

13. check local branch names:
git branch

14. create a new branch locally and switch to it:
git checkout -b main
