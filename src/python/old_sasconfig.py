import os

###     BEGIN SYSADMIN EDIT ###
###     BEGIN SYSADMIN EDIT ###
###     BEGIN SYSADMIN EDIT ###

#arch = "cluster"
arch = "mac"
#arch = "linux"

__cuda__ = True
__cuda__ = False

__openmm__ = False

__level__ = 'WARNING'
__level__ = 'DEBUG'

###     END SYSADMIN EDIT ###
###     END SYSADMIN EDIT ###
###     END SYSADMIN EDIT ###


if arch == "cluster":

    installation_bin_path = ['share','apps','local','bin']
    __core_libraries_include__ = [os.path.join(os.path.sep, 'share', 'apps', 'local', 'core_libraries', 'include')]
    __core_libraries_lib__ = [os.path.join(os.path.sep, 'share', 'apps', 'local', 'core_libraries', 'lib')]

    if __cuda__:
        __cuda_path__ = os.path.join(os.path.sep, 'share', 'apps', 'local', 'cuda-6.0')
else:

    installation_bin_path = ['usr','local','bin']
    __core_libraries_include__ = [os.path.join(os.path.sep, 'usr', 'local', 'core_libraries', 'include')]
    __core_libraries_lib__ = [os.path.join(os.path.sep, 'usr', 'local', 'core_libraries', 'lib')]

    if __cuda__:
        __cuda_path__ = os.path.join(os.path.sep, 'usr', 'local', 'cuda')

__bin_path__ = os.path.join(os.path.sep,*installation_bin_path)



