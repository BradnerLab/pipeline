import os
import errno
import logging
import sys

from os.path import dirname 
from os.path import basename

version = '1.4.0'

chromosome_name_length = 64 # Includes 1 for null terminator, so really max of 63 characters.
                            # Note that changing this value requires updating C++ code as well.

def mkdir_if_not_exists(directory):
    try:
        os.mkdir(directory)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def most_appropriate_executable_path(executable):
    # bamliquidator/motif_liquidator may be run by either a developer install from a git pipeline checkout,
    # or from a user install so that the exectuable is on the path.  If it seems like this is a git checkout,
    # then we return the developer executable path, else just use standard path

    if basename(dirname(dirname(os.path.realpath(__file__)))) == 'bamliquidator_internal':
        # look for developer executable location 
        executable_path = os.path.join(dirname(dirname(os.path.realpath(__file__))), executable)
        if not os.path.isfile(executable_path):
            exit("%s is missing -- try cd'ing into the directory and running 'make'" % executable_path)
        return executable_path
    else:
        # just look on standard path
        return executable

def configure_logging(args, output_directory, quiet):
    # Using root logger so we can just do logging.info/warn/error in this and other files.
    # If people start using bamliquidator_batch as an imported module, then we should probably
    # change this logging to not use the root logger directly.
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    file_handler = logging.FileHandler(os.path.join(output_directory, 'log.txt'))
    file_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s\t%(message)s',
                                                datefmt='%Y-%m-%d %H:%M:%S'))
    
    logger.addHandler(file_handler)
    # todo: add bamliquidator version to the starting log message
    logging.info("Starting %s %s with args %s", basename(sys.argv[0]), version, vars(args))

    # Adding console handler after writing the startup log entry.  The startup log could be useful 
    # in a file that is being appended to from a prior run, but would be annonying on stderr.

    console_handler = logging.StreamHandler()
    if quiet:
        console_handler.setLevel(logging.ERROR)
    else:
        console_handler.setLevel(logging.INFO)

    class FormatterNotFormattingInfo(logging.Formatter):
        def __init__(self, fmt):
            logging.Formatter.__init__(self, fmt)

        def format(self, record):
            if record.levelno == logging.INFO:
                return record.getMessage()
            return logging.Formatter.format(self, record)

    console_handler.setFormatter(FormatterNotFormattingInfo('%(levelname)s\t%(message)s'))
    logger.addHandler(console_handler)

