import feynman
import subprocess
import os, sys

def main():
    path = feynman.__path__[0]
    command = os.path.join(path, "_feynman_symbolic_regress1")
    subprocess.call([command] + sys.argv[1:])

if __name__ == '__main__':
    main()
