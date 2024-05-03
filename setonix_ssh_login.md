## ssh login to setonix using keys
https://pawsey.atlassian.net/wiki/spaces/US/pages/51925870/Use+of+SSH+Keys+for+Authentication
```bash
##Linux/macOS
ssh-keygen -t ecdsa -b 521 -f ~/.ssh/pawsey_ecdsa_key

##Windows
ssh-keygen -t ecdsa -b 521 -f $env:USERPROFILE/.ssh/pawsey_ecdsa_key

### you are prompted to choose and type a passphrase to protect the private key from being used by whoever gets access to it

Generating public/private ecdsa key pair.
Enter passphrase (empty for no passphrase): [Type a passphrase]
Enter same passphrase again: [Type the passphrase again]

Your identification has been saved in /home/<user>/.ssh/pawsey_ecdsa_key.
Your public key has been saved in /home/<user>/.ssh/pawsey_ecdsa_key.pub.
The key fingerprint is:
SHA256:K8R/F6+nBeDNpRskOfl/FnwpTPiWI3WBPpbeHTMU8uk dip008@apple-kf
The key's randomart image is:
+---[ECDSA 521]---+
|             ..o.|
|           o..o.o|
|          *.oo++ |
|     .   . O=B=+.|
|      o S ..@BoE*|
|     . . .  oOo.+|
|      . o . o + o|
|       . . . o.o |
|            oo   |
+----[SHA256]-----+

##The user now has a public key, the pawsey_ecdsa_key.pub file, and a private key, the pawsey_ecdsa_key  file.
##Linux/macOS
ls ~/.ssh
pawsey_ecdsa_key
pawsey_ecdsa_key.pub

##Windows
dir $env:USERPROFILE\.ssh
 
    Directory: C:\Users\[username]\.ssh
 
Mode              LastWriteTime         Length Name
----              -------------         ------ ----
-a----      15/04/2021  9:40 AM           1766 pawsey_ecdsa_key
-a----      15/04/2021  9:40 AM            402 pawsey_ecdsa_key.pub

###Adding the private key to the SSH agent
##Linux/macOS
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/pawsey_ecdsa_key

##Next, modify the ~/.ssh/config file to automatically load keys into the SSH agent and store passphrases in the Keychain
Host *
 AddKeysToAgent yes
 UseKeychain yes
 IdentityFile ~/.ssh/pawsey_ecdsa_key

###Copy your public key to the server
##Linux/macOS
ssh-copy-id -i ~/.ssh/pawsey_ecdsa_key.pub <username>@<remotehost>

After setting up the ssh-key and the ssh-agent, access to Pawsey systems will not ask for your password, but use the ssh-key instead:

```
