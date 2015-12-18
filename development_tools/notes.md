
### 0. Identation

Here is an example snippet from a .vimrc file that defines indentation rules. Note that the entire code base will have to be filtered to implement this. Any new code should follow these rules. Developers that use other IDEs should provide their settings for others to benefit.

filetype plugin on
set smartindent
set ts=4
set shiftwidth=4
set expandtab

Here is the indentation rules for Wingware Python IDE (which has also been followed in Emacs24):

Default Tab Size: 4
Default Indent Size: 4


### SETUP.PY ###

# format for data files:  ( path in site-packages, [ path + file_name of file to move] )
# format for data files:  ( path in site-packages, [ path1 + file_name1 of file to move, path2 + filename2 of file to move] )
# I will list single files even if multiple files go to the same place so that editing is easier

    data_files = [
            ( os.path.join('sasmol','test_sasmol','manual_tests') , [os.path.join('src','python','test_sasmol','manual_tests','hiv1_gag.pdb')]),
                    ( os.path.join('sasmol','test_sasmol','manual_tests') , [os.path.join('src','python','test_sasmol','manual_tests','hiv1_gag_200_frames.dcd')])
                                       ]



### GIT ###

https://git-scm.com/book/en/v2/Customizing-Git-Git-Configuration

# config:

git config --list
git config --global user.name "Joseph E. Curtis"
git config --global user.email "madscatt@gmail.com"
git config --global core.editor "vi"

(jc@hal)git_working_copies/sasmol% cat ~/.gitignore_global 
*~
.DS_Store
*.so
*.pyc
*.o

git config --global core.excludesfile ~/.gitignore_global 

# both of these do not seem to work ???

git config --global core.autocrlf = true
git config --global core.autocrlf = input

# but if you merely edit the ~/.gitconfig file it works
# (there are other things in this file ... only showing the autocrlf bit

[core]
    autocrlf = input

# cloning


git clone https://github.com/madscatt/sasmol.git


git checkout -b add_cpp_section

... do stuff ...

git add ...

git commit -m ''

git push -u sasmol add_cpp_section


(then you CAN use GitHub to merge and delete branch)

(but afterwards on local you must

# update the local files already done on GitHub
git pull

# make the master the active branch
git checkout master

# delete the other branch
git branch -d add_cpp_section

#list branches
git branch -av
     
