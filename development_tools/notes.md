
### SETUP.PY ###

# format for data files:  ( path in site-packages, [ path + file_name of file to move] )
# format for data files:  ( path in site-packages, [ path1 + file_name1 of file to move, path2 + filename2 of file to move] )
# I will list single files even if multiple files go to the same place so that editing is easier

    data_files = [
            ( os.path.join('sasmol','test_sasmol','manual_tests') , [os.path.join('src','python','test_sasmol','manual_tests','hiv1_gag.pdb')]),
                    ( os.path.join('sasmol','test_sasmol','manual_tests') , [os.path.join('src','python','test_sasmol','manual_tests','hiv1_gag_200_frames.dcd')])
                                       ]



### GIT ###

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
     
