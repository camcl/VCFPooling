import os
import pathlib
import shutil

from src.VCFPooling.python.archived.myssh import sshtools

from persotools.files import dir_size

'''
Copy the local directories tree and chosen files to the Uppmax server.
To be run from local, NOT from Uppmax server.
calling environment's ssh configuration settings (~/.ssh/config)
'''
# TODO: create function for single file transfer?

# Steps
create_rep = True
add_venv = False
add_file = False

# $HOME paths settings
cwd = os.getcwd()
local = pathlib.PurePath('/home/camille/1000Genomes')
littlelocal = pathlib.PurePath('/home/camille/PycharmProjects/1000Genomes/')

print('\n'.ljust(80, '*'))
print('Creating local light repository and transfer it to the remote server')
if create_rep:
    ext_out = shutil.ignore_patterns('*.vcf', 'IMP*.gz', 'REF*.gz', 'REF*.csi', 'IMP*.csi',
                                     '*.png', '*.pdf',
                                     '*python*', '*review*', '*simpool*', '*archive*',
                                     '*venv*', '*ppmax*', '.idea')
    try:
        tree = shutil.copytree(local,
                               littlelocal,
                               ignore=ext_out)
    except FileExistsError:
        pass
    size_in = dir_size(littlelocal)

    sshtools.ssh_recursive_transfer(sshtools.host,
                                    str(littlelocal),  # + '/src',
                                    str(sshtools.rackham),  # + '/1000Genomes',
                                    username=sshtools.user,
                                    password=sshtools.pwd)

    # with ssh_connect(host=sshtools.host, user=sshtools.user, connect_kwargs={'key_filename': sshtools.ssh_key,}:
    # with sshtools.ssh_connect(host=sshtools.host, user=sshtools.user, connect_kwargs={'password': sshtools.pwd}) as client:
    #     pytest = client.run(sshtools.py3cmd(['import os',
    #                                          'cwd = os.getcwd()',
    #                                          'print(cwd)']),
    #                         pty=True)

    sshtools.ssh_recursive_transfer(sshtools.host,
                                    littlelocal,
                                    sshtools.domus,
                                    username=sshtools.user,
                                    password=sshtools.pwd)
    shutil.rmtree(littlelocal)

    print('Bytes in --> ', size_in)
    print('Bytes out --> ', 'du -sh {}/1000Genomes'.format(sshtools.rackham))

print('\n'.ljust(80, '*'))
print('Adding venv: {}'.format(add_venv))
print(os.path.exists(pathlib.PurePath(sshtools.rackham, '/1000Genomes/venv3.6')))
if add_venv and os.path.exists(pathlib.PurePath(sshtools.rackham, '/1000G/venv3.6')):
    print('Add persotools library')
    sshtools.ssh_recursive_transfer(sshtools.host,
                                    '/home/camille/PycharmProjects/lib/persotools',
                                    str(sshtools.rackham) + '/1000Genomes/venv3.6/lib/python3.6/site-packages',
                                    username=sshtools.user,
                                    password=sshtools.pwd)