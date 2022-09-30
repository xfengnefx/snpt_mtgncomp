import gzip, sys
import binascii


def is_gz(filename):
    try:
        file = open(filename, 'rb')
    except:
        print('[warning::is_gz] File ({0}) not found.'.format(filename), file=sys.stderr)
        raise ValueError
    if binascii.hexlify(file.read(2)) == b'1f8b':
        file.close()
        return 0
    else:
        file.close()
        return 1


def opener(filename):
    status = is_gz(filename)
    if status==0: gzipped = True
    elif status==1: gzipped = False
    else: raise ValueError
    if gzipped:
        with gzip.open(filename, 'rb') as file:
            for line in file:
                yield line.decode()
    else:
        with open(filename) as file:
            for line in file:
                yield line


