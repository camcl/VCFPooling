import os
import pathlib
from contextlib import contextmanager
import paramiko
import fabric

from typing import *

from scripts.VCFPooling.poolSNPs.myssh import sshsession

'''
Tools for:
* establishing a SSH connection (Rackham server at Uppmax)
* create/load SSH keys
* transfer files with scp
* ...
'''

# $HOME paths settings
uppmax = 'camcl609@rackham.uppmax.uu.se'
domus = pathlib.PurePath('/domus/h1/camcl609')
rackham = pathlib.PurePath('/home/camcl609')
# generate key
ssh_key = pathlib.PurePath('/home/camille/.ssh/id_rackham')
# in terminal, hash the hosts file:
# $ ssh-keygen -H -f ~/.ssh/known_hosts
host_key_file = pathlib.PurePath('/home/camille/.ssh/known_hosts')

# SSH connection
host = 'rackham.uppmax.uu.se'
ip_host = '89.44.250.83'
user = 'camcl609'
pwd = 'Meew9iej31285'

hkey = paramiko.hostkeys.HostKeys(filename=str(host_key_file))
hkey.save(host_key_file)
# Obligatory warning: Do not use AutoAddPolicy – You are losing a protection against MITM attacks by doing so.
# For a correct solution, see Paramiko "Unknown Server".


@contextmanager
def ssh_connect(**kwargs):
    clt = fabric.Connection(**kwargs)
    yield clt
    clt.close()


def py3cmd(cmd_list: List[str]) -> str:
    pycmd = 'python3 -c' + ' "' + '; '.join(cmd_list) + '"'
    return pycmd


def get_cmd_output(clt: paramiko.SSHClient, command: str) -> tuple:
    stdin, stdout, stderr = clt.exec_command(command, get_pty=True)
    out = stdout.read()
    err = stderr.read()
    return out, err


# SSH/SCP Directory Recursively
def ssh_recursive_transfer(hostname,  src: Union[str, pathlib.Path], dst: Union[str, pathlib.Path],
                           username='root', key_file=None, password=None):
    session = sshsession.SSHSession(hostname, username=username, key_file=key_file, password=password)
    if not os.path.exists(dst):
        print('mkdir ' + str(dst))
        session.run('mkdir {}'.format(dst))
    session.put_all(src, dst)

