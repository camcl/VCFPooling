import os
import pathlib
import shutil

from scripts.poolSNPs.myssh import sshtools

from persotools.files import dir_size

'''
Copy the local directories tree and chosen files to the Uppmax server.
To be run from local, NOT from Uppmax server.
calling environment's ssh configuration settings (~/.ssh/config)
'''

# Steps
create_rep = False
add_venv = True
add_file = False

# $HOME paths settings
cwd = os.getcwd()
local = pathlib.PurePath('/home/camille/1000Genomes')
littlelocal = pathlib.PurePath('/home/camille/PycharmProjects/1000G/')

if create_rep:
    ext_out = shutil.ignore_patterns('*.vcf', 'IMP*.gz', 'REF*.gz', 'REF*.csi', 'IMP*.csi',
                                     '*.png', '*.pdf',
                                     '*ppmax*', '*review*', '*simpool*', '*archive*',
                                     '*venv*')
    try:
        tree = shutil.copytree(local,
                               littlelocal,
                               ignore=ext_out)
    except FileExistsError:
        pass
    size_in = dir_size(littlelocal)


    sshtools.ssh_recursive_transfer(sshtools.host,
                                    str(littlelocal) + '/scripts',
                                    str(sshtools.rackham) + '/1000G',
                                    username=sshtools.user,
                                    password=sshtools.pwd)

    # with ssh_connect(host=sshtools.host, user=sshtools.user, connect_kwargs={'key_filename': sshtools.ssh_key,}:
    with sshtools.ssh_connect(host=sshtools.host, user=sshtools.user, connect_kwargs={'password': sshtools.pwd}) as client:
        pytest = client.run(sshtools.py3cmd(['import os',
                                             'cwd = os.getcwd()',
                                             'print(cwd)']),
                            pty=True)
        # t = fabric.connection.Transfer(client)
        # t.put(str(local) + '/data/adaptative_gl.csv', remote=str(rackham))

    sshtools.ssh_recursive_transfer(sshtools.host,
                                    littlelocal,
                                    sshtools.domus,
                                    username=sshtools.user,
                                    password=sshtools.pwd)
    shutil.rmtree(littlelocal)

    with sshtools.ssh_connect(host=sshtools.host, user=sshtools.user, connect_kwargs={'password': sshtools.pwd}) as client:
        ls = client.run('ls', pty=True)
        size_out = client.run('du -sh {}/1000G'.format(sshtools.rackham))

    print('\n')
    print('Bytes in --> ', size_in)
    print('\n')
    print('Bytes out --> ', size_out)

print(add_venv)
print(os.path.exists(pathlib.PurePath(sshtools.rackham, '/1000G/venv')))
if add_venv and os.path.exists(pathlib.PurePath(sshtools.rackham, '/1000G/venv')):
    print('Add persotools library')
    sshtools.ssh_recursive_transfer(sshtools.host,
                                    '/home/camille/PycharmProjects/lib/persotools',
                                    str(sshtools.rackham) + '/1000G/venv/lib/python3.6/site-packages',
                                    username=sshtools.user,
                                    password=sshtools.pwd)
    sshtools.ssh_recursive_transfer(sshtools.host,
                                    '/home/camille/1000Genomes/scripts/__init__.py',
                                    str(sshtools.rackham) + '/1000G/venv/lib/python3.6/site-packages/scripts',
                                    username=sshtools.user,
                                    password=sshtools.pwd)
