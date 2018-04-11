from __future__ import print_function, division, absolute_import #, unicode_literals

import sys
import os
import time
import shutil
import textwrap
import platform
import tarfile
import re
import warnings

from socket import gethostname
from subprocess import Popen, PIPE

# Handle py2, py3k differences.
py2 = sys.version_info[0] <= 2
if py2:
    import cPickle as pickle
    from StringIO import StringIO
    from ConfigParser import SafeConfigParser, RawConfigParser, ConfigParser, NoOptionError
else:
    import pickle
    from io import StringIO
    from configparser import SafeConfigParser, RawConfigParser, ConfigParser, NoOptionError

from .jobrunner import TimeBomb
from .tools import (RestrictedShell,  StringColorizer, unzip, tail_file,
                    pprint_table, Patcher, Editor)
from .xyaptu import xcopier
from .devtools import FileLock
from .memprof import AbimemParser
from .termcolor import cprint

from collections import namedtuple
# OrderedDict was added in 2.7. ibm6 still uses python6
try:
    from collections import OrderedDict
except ImportError:
    from .ordereddict import OrderedDict

import logging
logger = logging.getLogger(__name__)

__version__ = "0.5"
__author__ = "Matteo Giantomassi"

__all__ = [
    "BuildEnvironment",
    "AbinitTestSuite",
]

_MY_NAME = os.path.basename(__file__)[:-3] + "-" + __version__


# Helper functions and tools

def fix_punctuation_marks(s):
    """
    Remove whitespaces before `,` and `;` that trigger a bug in ConfigParser
    when the option spans multiple lines separated by `;`

    For instance:

        #%% files_to_test =
        #%%   t93.out, tolnlines = 0, tolabs = 0.000e+00, tolrel = 0.000e+00 ;
        #%%   t93.out_ep_SBK, tolnlines = 0, tolabs = 0.000e+00, tolrel = 0.000e+00

    is not treated properly by ConfigParser due to ` ;` at the end of the line and
    only t93.out is stored in dictionary and tested by fldiff.
    """
    for mark in (",", ";"):
        s = (mark + " ").join([tok.rstrip() for tok in s.split(mark)])
    return s + "\n"


def html_colorize_text(string, code):
    return "<FONT COLOR='%s'>%s</FONT>" % (code, string)

_status2htmlcolor = {
    "succeeded": lambda string : html_colorize_text(string, 'Green'),
    "passed": lambda string : html_colorize_text(string, 'DeepSkyBlue'),
    "failed": lambda string : html_colorize_text(string, 'Red'),
    "disabled": lambda string : html_colorize_text(string, 'Cyan'),
    "skipped": lambda string : html_colorize_text(string, 'Cyan'),
}


def status2html(status):
    """Convert test status in a colored HTML string."""
    return _status2htmlcolor[status](status)


def sec2str(seconds):
    """Convert seconds to string."""
    return "%.2f" % seconds


def str2html(string, end="<br>"):
    """Returns a HTML string."""
    lines = string.splitlines()
    return "<br>".join(lines) + end


def args2htmltr(*args):
    string = ""
    for arg in args:
        string += "<td>" + str(arg) + "</td>"
    return string


def html_link(string, href=None):
    """Create a HTML link from a string. Use href as link of href is not None."""
    if href is not None:
        return "<a href='%s'>%s</a>" % (href, string)
    else:
        return "<a href='%s'>%s</a>" % (string, string)


def is_string(s):
    try:
        s + "hello"
        return True
    except TypeError:
        return False


def has_exts(path, exts):
    """True if path ends with extensions exts"""
    root, ext = os.path.splitext(path)
    if is_string(exts):
        return ext == exts
    else:
        return ext in exts


def lazy__str__(func):
    """Lazy decorator for __str__ methods"""
    def oncall(*args, **kwargs):
        self = args[0]
        return "\n".join([str(k) + " : " + str(v) for (k, v) in self.__dict__.items()])
    return oncall

# Helper functions for performing IO


def lazy_read(fname):
    if sys.version_info >= (3, 0):
        #with open(fname, "rt", encoding="ISO-8859-1") as fh:
        #with open(fname, "rt", encoding="utf-8", errors="ignore") as fh:
        with open(fname, "rt", encoding="utf-8") as fh:
            return fh.read()

    else:
        with open(fname, "rt") as fh:
            return fh.read()


def lazy_readlines(fname):
    if sys.version_info >= (3, 0):
        with open(fname, "rt", encoding="utf-8") as fh:
            return fh.readlines()
    else:
        with open(fname, "rt") as fh:
            return fh.readlines()


def lazy_write(fname, s):
    with open(fname, "wt") as fh:
        fh.write(s)


def lazy_writelines(fname, lines):
    with open(fname, "wt") as fh:
        fh.writelines(lines)


class Record(object):
    @lazy__str__
    def __str__(self):
        pass


def rmrf(top, exclude_paths=None):
    """
    Recursively remove all files and directories contained in directory top.

    Args:
        exclude_paths:
            list with the absolute paths that should be preserved

    Returns the list of files and the directories that have been removed.
    """
    exc_paths = []
    if exclude_paths is not None:
        if is_string(exclude_paths):
            exc_paths = [exclude_paths]
        else:
            exc_paths = exclude_paths

    removed = []
    for (root, dirs, files) in os.walk(top):
        for f in files:
            file_path = os.path.join(root, f)
            if file_path not in exc_paths:
                os.unlink(file_path)
                removed.append(file_path)
        for d in dirs:
            dir_path = os.path.join(root, d)
            if dir_path not in exc_paths:
                shutil.rmtree(dir_path)
                removed.append(dir_path)

    return removed


def find_abortfile(workdir):
    """
    Return the absolute path of the MPIABORTFILE file produced by (Abinit|Abinit_with_libpaw)
    Empty string if file is not present.

    Args:
        workdir: Working directory of the test.

    .. Note::

        __LIBPAW_MPIABORFILE__ is produced if abinit uses libpaw and when we die inside libpaw.
    """
    for s in ("__ABI_MPIABORTFILE__", "__LIBPAW_MPIABORFILE__"):
        path = os.path.join(workdir, s)
        if os.path.exists(s): return path
    return ""


def read_yaml_errmsg(path):
    """
    Extract the YAML error message from file `path`.
    Returns string with message, empty string if message is not found.

    The Yaml error message is in the form:

    --- !ERROR
    src_file: m_io_screening.F90
    src_line: 648
    message: |
        Unsupported value of iomode
    ...
    """
    errlines, inerr = [], 0

    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("---") and ("ERROR" in line or "BUG" in line):
                inerr = 1

            if inerr:
                errlines.append(line)
                if line.startswith("..."): break

    return "".join(errlines)


def extract_errinfo_from_files(workdir):
    """
    Extract information from the files produced by the code when we run tests in debug mode.

    Return:
        String with the content of the files. Empty string if no debug file is found.
    """
    registered_exts = set([".flun", ".mocc"])
    errinfo = []

    for path in os.listdir(workdir):
        _, ext = os.path.splitext(path)
        if ext not in registered_exts: continue
        with open(os.path.join(workdir, path), "rt") as fh:
            errinfo.append(" ")
            errinfo.append("From file: %s " % path)
            errinfo.extend(l.strip() for l in fh)

    return "\n".join(errinfo)


class FileToTest(object):
    """This object contains information on the output file that will be analyzed by fldiff"""
    #  atr_name,   default, conversion function. None designes mandatory attributes.
    _attrbs = [
        ("name",     None, str),
        ("tolnlines",None, int),    # fldiff tolerances
        ("tolabs",   None, float),
        ("tolrel",   None, float),
        ("fld_options","",str) ,     # options passed to fldiff.
        ("fldiff_fname","",str),
        ("hdiff_fname","",str),
        ("diff_fname","",str),
        #("pydiff_fname","",str),
    ]

    def __init__(self, dic):

        for atr in FileToTest._attrbs:
            atr_name = atr[0]
            default = atr[1]
            f = atr[2]
            value = dic.get(atr_name, default)
            if value is None:
                raise ValueError("%s must be defined" % atr_name)

            value = f(value)
            if hasattr(value, "strip"): value = value.strip()
            self.__dict__[atr_name] = value

        # Postprocess fld_options
        self.fld_options = self.fld_options.split()
        for opt in self.fld_options:
            if not opt.startswith("-"):
                raise ValueError("Wrong fldiff option: %s" % opt)

    @lazy__str__
    def __str__(self): pass

    def compare(self, fldiff_path, ref_dir, workdir, timebomb=None, outf=sys.stdout):
        """
        Use fldiff_path to compare the reference file located in ref_dir with
        the output file located in workdir. Results are written to stream outf.
        """
        ref_fname = os.path.abspath(os.path.join(ref_dir, self.name))
        # FIXME Hack due to the stdout-out ambiguity
        if not os.path.exists(ref_fname) and ref_fname.endswith(".stdout"):
            ref_fname = ref_fname[:-7] + ".out"
        out_fname = os.path.abspath(os.path.join(workdir, self.name))

        opts = self.fld_options
        label = self.name

        fld_result, got_summary = wrap_fldiff(fldiff_path, ref_fname, out_fname,
                                              opts=opts, label=label, timebomb=timebomb, out_filobj=outf)

        if not got_summary:
            # Wait 10 sec, then try again (workaround for woopy)
            logger.critical("Didn't got fldiff summary, will sleep for 10 s...")
            time.sleep(10)
            fld_result, got_summary = wrap_fldiff(fldiff_path, ref_fname, out_fname,
                                                  opts=opts, label=label, timebomb=timebomb, out_filobj=outf)

            if not got_summary:
                logger.critical("fldiff summary is still empty!")

        isok, status, msg = fld_result.passed_within_tols(self.tolnlines, self.tolabs, self.tolrel)

        # Save comparison results.
        self.fld_isok = isok
        self.fld_status = status
        self.fld_msg = msg

        return isok, status, msg

# Parsers used for the different TEST_INFO options


def _str2filestotest(string):
    """
    Parse the files_to_test section.
    Returns a tuple of `FileToTest` objects.
    """
    if not string:
        return []

    if ";" in string:
        file_lines = [s for s in string.split(";") if s.strip()]
    else:
        file_lines = [string]

    files_to_test = []
    for line in file_lines:
        tokens = line.split(",")
        d = {"name": tokens[0]}
        for tok in tokens[1:]:
            k, v = [s.strip() for s in tok.split("=")]
            if k in d:
                err_msg = "Found multiple occurences of keyword %s" % k
                raise AbinitTestInfoParserError(err_msg)
            d[k] = v
        files_to_test.append(FileToTest(d))

    return tuple(files_to_test)


def _str2list(string):    return [s.strip() for s in string.split(",") if s]
def _str2intlist(string): return [int(item) for item in _str2list(string) ]
def _str2set(string):     return set([s.strip() for s in string.split(",") if s])
def _str2cmds(string):    return [s.strip() for s in string.split(";") if s]

def _str2bool(string):
    string = string.strip().lower()
    if string == "yes": return True
    return False

# TEST_INFO specifications
TESTCNF_KEYWORDS = {
# keyword        : (parser, default, section, description)
# [setup]
"executable"     : (str       , None , "setup", "Name of the executable e.g. abinit"),
"test_chain"     : (_str2list , ""   , "setup", "Defines a ChainOfTest i.e. a list of tests that are connected together."),
"need_cpp_vars"  : (_str2set  , ""   , "setup", "CPP variables that must be defined in config.h in order to enable the test."),
"exclude_hosts"  : (_str2list , ""   , "setup", "The test is not executed if we are running on a slave that matches compiler@hostname"),
"exclude_builders": (_str2list, ""   , "setup", "The test is not executed if we are using a builder whose name is in the list"),
"input_prefix"   : (str       , ""   , "setup", "Prefix for input files (used for the ABINIT files file)"),
"output_prefix"  : (str       , ""   , "setup", "Prefix for output files (used for the ABINIT files file)"),
"expected_failure": (_str2bool, "no" , "setup", "yes if the subprocess executing executable is expected to fail (retcode != 0) (default: no)"),
"input_ddb"      : (str       , ""   , "setup", "The input DDB file read by anaddb"),
"input_gkk"      : (str       , ""   , "setup", "The input GKK file read by anaddb"),
"system_xml"     : (str       , ""   , "setup","The system.xml file read by multibinit"),
"coeff_xml"      : (str       , ""   , "setup","The coeff.xml file read by multibinit"),
"md_hist"        : (str       , ""   , "setup","The hist file file read by multibinit"),
# [files]
"files_to_test"  : (_str2filestotest, "", "files", "List with the output files that are be compared with the reference results. Format:\n" +
                                                   "\t file_name, tolnlines = int, tolabs = float, tolrel = float [,fld_options = -medium]\n" +
                                                   "\t    tolnlines: the tolerance on the number of differing lines\n" +
                                                   "\t    tolabs:the tolerance on the absolute error\n" +
                                                   "\t    tolrel: tolerance on the relative error\n" +
                                                   "\t    fld_options: options passed to fldiff.pl (optional).\n" +
                                                   "\t    Multiple files are separated by ; e.g.\n" +
                                                   "\t    foo.out, tolnlines = 2, tolabs = 0.1, tolrel = 1.0e-01;\n" +
                                                   "\t    bar.out, tolnlines = 4, tolabs = 0.0, tolrel = 1.0e-01"
                                                       ),
"psp_files"      : (_str2list,        "", "files", "List of pseudopotential files (located in the Psps_for_tests directory)."),
"extra_inputs"   : (_str2list,        "", "files", "List of extra input files."),
# [shell]
"pre_commands"   : (_str2cmds, "", "shell", "List of commands to execute before starting the test"),
"post_commands"  : (_str2cmds, "", "shell", "List of commands to execute after the test is completed"),
# [paral_info]
"max_nprocs"     : (int      ,  1 , "paral_info", "Maximum number of MPI processors (1 for sequential run)"),
"nprocs_to_test" : (_str2intlist, "","paral_info","List with the number of MPI processes that should be used for the test"),
"exclude_nprocs" : (_str2intlist, "","paral_info","List with the number of MPI processes that should not be used for the test"),
# [extra_info]
"authors"         : (_str2set , "Unknown"                 , "extra_info", "Author(s) of the test"),
"keywords"       : (_str2set , ""                         , "extra_info", "List of keywords associated to the test"),
"description"    : (str      , "No description available",  "extra_info", "String containing extra information on the test"),
"topics"         : (_str2list, "",  "extra_info", "Topics associated to the test"),
"references"     : (_str2list, "",  "extra_info", "List of references to papers or other articles"),
}

#TESTCNF_SECTIONS = set( [ TESTCNF_KEYWORDS[k][2] for k in TESTCNF_KEYWORDS ] )

# This extra list is hardcoded in order to have a fixed order of the sections in doc_testcfn_format.
# OrderedDict have been introduced in python2.7 sigh!
TESTCNF_SECTIONS = [
  "setup",
  "files",
  "shell",
  "paral_info",
  "extra_info",
]

# consistency check.
for key, tup in TESTCNF_KEYWORDS.items():
    if tup[2] not in TESTCNF_SECTIONS:
        raise ValueError("Please add the new section %s to TESTCNF_SECTIONS" % tup[2])


def line_starts_with_section_or_option(string):
    """True if string start with a TEST_INFO section or option."""
    from re import compile
    re_ncpu = compile("^NCPU_(\d+)$")
    s = string.strip()
    idx = s.find("=")
    if idx == -1: # might be a section.
        if s.startswith("[") and s.endswith("]"):
            if s[1:-1] in TESTCNF_SECTIONS: return 1 # [files]...
            if re_ncpu.search(s[1:-1]): return 1    # [NCPU_1] ...
    else:
        if s[:idx].strip() in TESTCNF_KEYWORDS: return 2

    return 0


def doc_testcnf_format(fh=sys.stdout):
    """Automatic documentation of the TEST_INFO sections and related options."""
    def writen(string): fh.write(string + "\n")

    writen("Automatic documentation of the TEST_INFO sections and options.")

    for section in TESTCNF_SECTIONS:
        writen("\n["+section+"]")
        for key in TESTCNF_KEYWORDS:
            tup = TESTCNF_KEYWORDS[key]
            if section == tup[2]:
                line_parser = tup[0]
                default = tup[1]
                if default is None:
                    default = "Mandatory"
                desc = tup[3]
                if default:
                    msg = "%s =  %s (DEFAULT: %s)" % (key, desc, default)
                else:
                    msg = "%s =  %s" % (key, desc)
                writen(msg)


class AbinitTestInfo(object):
    """Container storing the options specified in the TEST_INFO section."""
    def __init__(self, dct):
        for k, v in dct.items():
            self.__dict__[k] = v

        #if self.nprocs_to_test and self.test_chain:
        #  err_msg = "test_chain and nprocs_to_test are mutually exclusive"
        #  raise TestInfoParserError(err_msg)

        # Add the executable name to the list of keywords.
        self.add_keywords([self.executable])

    @lazy__str__
    def __str__(self): pass

    def add_cpp_vars(self, need_cpp_vars):
        """Add new set of CPP variables."""
        self.need_cpp_vars = self.need_cpp_vars.union(need_cpp_vars)

    def add_keywords(self, keywords):
        """Add new set of keywords."""
        self.keywords = self.keywords.union(keywords)

    def make_test_id(self):
        """
        Generate the string with the test identifier
        A special treatment is used for the multi-parallel tests.
        In this case, the test_id is constructed by appending the string _MPI#
        where # is the number of MPI processors.
        """
        ## FIXME Assumes inp_fname is in the form name.in
        test_id = os.path.basename(self.inp_fname).split(".")[0]
        if self.ismulti_parallel:
            test_id += "_MPI%d" % self.max_nprocs
        return test_id

    @property
    def ismulti_parallel(self):
        """True is this is a multi-parallel test."""
        return self._ismulti_paral


class AbinitTestInfoParserError(Exception):
    """Error class raised by the parse"""


class AbinitTestInfoParser(object):
    """This object parses the TEST_INFO section that describes the test."""
    Error = AbinitTestInfoParserError

    def __init__(self, inp_fname, defaults=None):
        """
        Args:
            inp_fname:
                test input file
            defaults:
                default values passed to the INI parser.
        """
        logger.info("Parsing TEST_INFO section from input file : " + str(inp_fname))

        self.inp_fname = os.path.abspath(inp_fname)
        self.inp_dir, x = os.path.split(self.inp_fname)

        SENTINEL = '#%%'
        HEADER = "<BEGIN TEST_INFO>\n"
        FOOTER = "<END TEST_INFO>\n"

        # Extract the lines that start with SENTINEL and parse the file.
        lines = lazy_readlines(inp_fname)
        #for l in lines: print(l)
        #print(inp_fname)
        lines = [l.replace(SENTINEL, "", 1).lstrip() for l in lines if l.startswith(SENTINEL)]

        try:
            start, stop = lines.index(HEADER), lines.index(FOOTER)
        except ValueError:
            raise self.Error("%s does not contain any valid testcnf section!" % inp_fname)

        lines = lines[start+1:stop]
        if not lines:
            raise self.Error("%s does not contain any valid testcnf section!" % inp_fname)

        # Hack to allow options occupying more than one line.
        string = ""
        for l in lines:
            # MGDEBUG
            # This is needed to avoid problems with multiple lines, see docstring.
            l = fix_punctuation_marks(l)
            if line_starts_with_section_or_option(l):
                string += l
            else:
                if l.startswith("#"): continue
                string = string.rstrip() + " " + l
        lines = [l + "\n" for l in string.split("\n")]
        #MGDEBUG
        #print("in gmatteo's parser ")
        #for l in lines: print(l, end="")

        s = StringIO()
        s.writelines(lines)
        s.seek(0)

        class MySafeConfigParser(SafeConfigParser):
            """Wrap the get method of SafeConfigParser to disable the interpolation of raw_options."""
            raw_options = ["description",]

            def get(self, section, option, raw=False, vars=None):
                if option in self.raw_options and section == TESTCNF_KEYWORDS[option][2]:
                    logger.debug("Disabling interpolation for section = %s, option = %s" % (section, option))
                    #print("Disabling interpolation for section = %s, option = %s" % (section, option))
                    if py2:
                        return SafeConfigParser.get(self, section, option, raw=True, vars=vars)
                    else:
                        return SafeConfigParser.get(self, section, option, raw=True, vars=vars, fallback=None)
                else:
                    #print("Calling SafeConfigParser for section = %s, option = %s" % (section, option))
                    if py2:
                        return SafeConfigParser.get(self, section, option, raw, vars)
                    else:
                        return SafeConfigParser.get(self, section, option, raw=raw, vars=vars, fallback=None)

        # Old version. ok with py2 but not with py3k
        if py2:
            self.parser = MySafeConfigParser(defaults) # Wrap the parser.
            #self.parser = RawConfigParser(defaults)
        else:
            self.parser = ConfigParser(defaults, interpolation=None)

        try:
            self.parser.readfp(s)
        except Exception as exc:
            cprint("Exception while parsing: %s\n%s" % (inp_fname, exc), "red")
            for l in lines: print(l, end="")
            raise exc

        # Consistency check
        opt = "test_chain"
        section = TESTCNF_KEYWORDS[opt][2]
        pars = TESTCNF_KEYWORDS[opt][0]

        if self.parser.has_option(section, opt):
            string = self.parser.get(section, opt)
            chain = pars(string)
            ones = [chain.count(value) for value in chain]
            if sum(ones) != len(ones):
                err_msg = "%s : test_chain contains repeated tests %s" % (inp_fname, string)
                raise self.Error(err_msg)

            # Check whether (section, option) is correct.
            #_defs = [s.upper() for s in defaults] if defaults else []
            #err_msg = ""
            #for section in parser.sections():
            #  for opt in parser.options(section):
            #    if opt.upper() in _defs: continue
            #    if opt not in TESTCNF_KEYWORDS:
            #      err_msg += "Unknown (section, option) = %s, %s\n" % (section, opt)
            #    elif section != TESTCNF_KEYWORDS[opt][2]:
            #      err_msg += "Wrong (section, option) = %s, %s\n" % (section, opt)
            #if err_msg: raise ValueError(err_msg)

    def generate_testinfo_nprocs(self, nprocs):
        """Returns a record with the variables needed to handle the job with nprocs."""
        info = Record()
        d = info.__dict__

        # First read and parse the global options.
        for key in TESTCNF_KEYWORDS:
            tup = TESTCNF_KEYWORDS[key]
            line_parser = tup[0]
            section = tup[2]

            if section in self.parser.sections():
                try:
                    d[key] = self.parser.get(section, key)
                except NoOptionError:
                    d[key] = tup[1] # Section exists but option is not specified. Use default value.
            else:
                d[key] = tup[1] # Section does not exist. Use default value.

            # Process the line
            #if key == "files_to_test": print("hello files", d[key], "bye files")
            try:
                d[key] = line_parser(d[key])
            except Exception as exc:
                try:
                    err_msg = "Wrong line:\n key = %s, d[key] = %s\n in file: %s" % (key, d[key], self.inp_fname)
                except:
                    err_msg = "In file %s:\n%s" % (self.inp_fname, str(exc))

                raise self.Error(err_msg)

        # At this point info contains the parsed global values.
        # Now check if this is a parallel test and, in case, overwrite the values
        # using those reported in the [CPU_nprocs] sections.
        # Set also the value of info._ismulti_paral so that we know how to create the test id
        if not info.nprocs_to_test:
            assert nprocs == 1
            info._ismulti_paral = False
        else:
            logger.debug("multi parallel case")
            if nprocs not in info.nprocs_to_test:
                err_msg = "in file: %s. nprocs = %s > not in nprocs_to_test = %s" % (self.inp_fname, nprocs, info.nprocs_to_test)
                raise self.Error(err_msg)

            if nprocs > info.max_nprocs:
                try:
                    err_msg = "in file: %s. nprocs = %s > max_nprocs = %s" % (self.inp_fname, nprocs, self.max_nprocs)
                except Exception as exc:
                    err_msg = "in file: %s\n%s" % (self.inp_fname, str(exc))

                raise self.Error(err_msg)

            # Redefine variables related to the number of CPUs.
            info._ismulti_paral = True
            info.nprocs_to_test = [nprocs]
            info.max_nprocs = nprocs

            info.exclude_nprocs = list(range(1, nprocs))
            #print(self.inp_fname, nprocs, info.exclude_nprocs)

            ncpu_section = "NCPU_" + str(nprocs)
            if not self.parser.has_section(ncpu_section):
                err_msg = "Cannot find section %s in %s" % (ncpu_section, self.inp_fname)
                raise self.Error(err_msg)

            for key in self.parser.options(ncpu_section):
                if key in self.parser.defaults(): continue
                opt = self.parser.get(ncpu_section, key)
                tup = TESTCNF_KEYWORDS[key]
                line_parser = tup[0]
                #
                # Process the line and replace the global value.
                try:
                    d[key] = line_parser(opt)
                except:
                    err_msg = "In file: %s. Wrong line: key: %s, value: %s" % (self.inp_fname, key, d[key])
                    raise self.Error(err_msg)

                #print(self.inp_fname, d["max_nprocs"])

        # Add the name of the input file.
        info.inp_fname = self.inp_fname

        return AbinitTestInfo(d)

    @property
    def nprocs_to_test(self):
        """List with the number of MPI processors to be tested."""
        key = "nprocs_to_test"
        opt_parser = TESTCNF_KEYWORDS[key][0]
        default = TESTCNF_KEYWORDS[key][1]
        section = TESTCNF_KEYWORDS[key][2]

        try:
            opt = self.parser.get(section, key)
        except NoOptionError:
            opt = default

        return opt_parser(opt)

    @property
    def is_testchain(self):
        """True if this is a chain of tests"""
        opt = "test_chain"
        section = TESTCNF_KEYWORDS[opt][2]
        return self.parser.has_option(section, opt)

    def chain_inputs(self):
        """Return a list with the path of the input files belonging to the test chain"""
        assert self.is_testchain
        opt = "test_chain"
        section = TESTCNF_KEYWORDS[opt][2]
        parse = TESTCNF_KEYWORDS[opt][0]

        fnames = parse(self.parser.get(section, opt))
        return [os.path.join(self.inp_dir, fname) for fname in fnames]

    #@property
    #def is_parametrized_test(self):
    #    """True if this is a parametrized test."""
    #    raise NotImplemented()

    #def get_parametrized_tests(self):
    #    """Return the list of parametrized tests."""


def find_top_build_tree(start_path, with_abinit=True, ntrials=10):
    """
    Returns the absolute path of the ABINIT build tree.
    Assume start_path is within the build tree.

    Raises:
        `RuntimeError` if build tree is not found after ntrials attempts.
    """
    abs_path = os.path.abspath(start_path)

    trial = 0
    while trial <= ntrials:
        config_h = os.path.join(abs_path, "config.h")
        abinit_bin = os.path.join(abs_path, "src", "98_main", "abinit")
        # Check if we are in the top of the ABINIT source tree
        if with_abinit:
            found = os.path.isfile(config_h) and os.path.isfile(abinit_bin)
        else:
            found = os.path.isfile(config_h)

        if found:
            return abs_path
        else:
            abs_path, tail = os.path.split(abs_path)
            trial += 1

    raise RuntimeError("Cannot find the ABINIT build tree after %s trials" % ntrials)


class Compiler(object):
    """
    Base class for C,Fortran,C++ compilers.
    Usually instantiated through the class method from_defined_cpp_vars.
    """
    def __init__(self, name, version=None):
        self.name = name
        self.version = version

    def __str__(self):
        return "%s: %s %s" % (self.__class__.__name__, self.name, self.version)

    @classmethod
    def from_defined_cpp_vars(cls, defined_cpp_vars):
        for var in defined_cpp_vars:
            # TODO: version may be useful but it's not reported in config.h
            if var in cls._KNOWN_CPP_VARS:
                # Build the name of the compiler.
                name = var.lower().split("_")[1]
                if name == "gnu": name = "gfortran"
                if name == "pathscale": name = "psc"
                return cls(name=name, version=None)
        else:
            err_msg = "Cannot detect the name of the %s\n. Defined CPP vars: %s " % (cls.__name__, str(defined_cpp_vars))
            raise RuntimeError(err_msg)


class FortranCompiler(Compiler):
    """
    Store information on the Fortran compiler used to build abinit.
    """
    # CPP variables used in config.h
    _KNOWN_CPP_VARS = [
        "FC_ABSOFT",
        "FC_FUJITSU",
        "FC_G95",
        "FC_GNU",
        "FC_HITACHI",
        "FC_IBM",
        "FC_INTEL",
        "FC_MIPSPRO",
        "FC_NAG",
        "FC_OPEN64",
        "FC_PATHSCALE",
        "FC_PGI",
        "FC_SUN",
    ]


class CPreProcessorError(Exception):
    """Errors raised by `CPreProcessors`"""


class CPreProcessor(object):
    """Pre-process source code with ANSI CPP."""
    Error = CPreProcessorError

    def __init__(self, includes=None, opts=None, bin="cpp", verbose=0):
        self.includes = ["."]
        if includes is not None: self.includes = includes
        self.opts = ["-DHAVE_CONFIG_H"]
        if opts is not None: self.opts = opts
        self.bin, self.verbose = bin, verbose

    def process_file(self, filepath, remove_lhash=True):
        """
        Read source from filepath, call CPP wit the includes and the
        options passed to the constructor.

        Returns:
            preprocessed text.
        """
        if self.bin is None:
            # No pre-processing, return raw string.
            with open(filepath, "r") as f:
                return f.read()

        cmd = [self.bin]
        if self.opts: cmd += self.opts
        cmd += ["-ansi"]
        if self.includes: cmd += ["-I"+inc for inc in self.includes]
        cmd += [filepath]
        cmd = " ".join(cmd)
        if self.verbose: print(cmd)

        from subprocess import Popen, PIPE
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        if p.returncode:
            raise self.Error("C-preprocessor returned %d\n stderr:\n%s" % (p.returncode, stderr))

        # Remove leading hash symbols added by CPP
        if not remove_lhash:
            return stdout
        else:
            return "\n".join([str(l) for l in stdout.splitlines() if not l.startswith("#")])


class FortranBacktrace(object):
    def __init__(self, text):
        self.text = text
        self.trace = []
        self.parse()

    def __str__(self):
        return str(self.trace)

    def parse(self):
        raise NotImplementedError("parse method must be implemented by the subclass")

    def locate_srcfile(self, base_name):
        top = find_top_build_tree(start_path=".", with_abinit=True)
        top = os.path.join(top, "src")

        for dirpath, dirnames, filenames in os.walk(top):
            if base_name in filenames:
                apath = os.path.join(dirpath, base_name)
                return apath
        else:
            cprint("Cannot find file: %s" % base_name, "red")
            return None

    def edit_source(self, editor=None):
        if not self.trace: return

        if editor is None: editor = Editor()
        src_file, lineno = self.trace[0]
        src_file = self.locate_srcfile(src_file)

        return editor.edit_file(src_file, lineno=lineno)


class NagBacktrace(FortranBacktrace):

    def parse(self):
        # Example
        #
        # Runtime Error: opernl4a_cpp.f90, line 871: INTEGER(int32) overflow for 2146435072 * 3
        # Program terminated by fatal error
        # opernl4a_cpp.f90, line 871: Error occurred in OPERNL4A
        if not self.text: return

        #MAGIC = "Program terminated by fatal error"
        #for i, line in enumerate(self.text):
        #   if MAGIC in line: break
        #else:
        #    return

        re_nagline = re.compile("(\w+\.f90), line (\d+): (.+)")

        for line in self.text:
            m = re_nagline.match(line)
            if not m: continue
            src_file, lineno = m.group(1), m.group(2)
            self.trace.append((src_file, int(lineno)))


class BuildEnvironment(object):
    """Store information on the build environment."""

    def __init__(self, build_dir, cygwin_instdir=None):
        """
        Args:
            build_dir:
                Path to the top level directory of the build.
            cygwin_instdir:
                Installation directory of cygwin. Defaults to '/cygwin'
        """
        # Try to figure out the top level directory of the build tree.
        try:
            build_dir = find_top_build_tree(build_dir)
        except:
            raise

        self.uname = platform.uname()
        self.hostname = gethostname().split(".")[0]

        try:
            self.username = os.getlogin()
        except:
            self.username = "No_username"

        self.build_dir = os.path.abspath(build_dir)
        self.configh_path = os.path.join(self.build_dir, "config.h")
        self.binary_dir = os.path.join(self.build_dir, "src", "98_main")

        self._cygwin_instdir = ""
        if cygwin_instdir is not None:
            self._cygwin_instdir = cygwin_instdir

        # Binaries that are not located in src/98_main
        self._external_bins = {
            "atompaw": os.path.join(self.build_dir, "fallbacks", "exports", "bin", "atompaw-abinit"),
            "timeout": os.path.join(self.build_dir, "tests", "Timeout", "timeout"),
        }

        # Check if this is a valid ABINIT build tree.
        if not (os.path.isfile(self.configh_path) and os.path.isfile(self.path_of_bin("abinit"))):
            raise ValueError("%s is not a valid ABINIT build tree." % self.build_dir)

        # Get the list of CPP variables defined in the build.
        self.defined_cppvars = parse_configh_file(self.configh_path)

        # Get info on the compilers
        self.fortran_compiler = FortranCompiler.from_defined_cpp_vars(self.defined_cppvars)
        #print(self.fortran_compiler)
        #if not self.has_bin("timeout"): print("Cannot find timeout executable!")

        self.buildbot_builder = None

    @lazy__str__
    def __str__(self): pass

    def issrctree(self):
        """True if this is a source tree."""
        configac_path = os.path.join(self.build_dir, "configure.ac")
        abinitF90_path = os.path.join(self.build_dir, "src", "98_main", "abinit.F90")

        return os.path.isfile(configac_path) and os.path.isfile(abinitF90_path)

    def iscygwin(self):
        """True if we are running under CYGWIN"""
        return "CYGWIN" in self.uname[0].upper()

    def _addext(self, string):
        """Append .exe extension, needed for cygwin"""
        if self.iscygwin(): string += ".exe"
        return string

    def path_of_bin(self, bin_name, try_syspath=True):
        """Return the absolute path of bin_name."""
        if bin_name in self._external_bins:
            bin_path = self._external_bins[bin_name]
        else:
            bin_path = os.path.join(self.binary_dir, bin_name) # It's in src/98_main

        bin_path = self._addext(bin_path)

        # Handle external bins that are installed system wide (such as atompaw on woopy)
        if bin_name in self._external_bins and not os.path.isfile(bin_path):
            if not try_syspath: return ""
            # Search it in PATH.
            paths = os.getenv("PATH").split(os.pathsep)
            for p in paths:
                bin_path = os.path.join(p, bin_name)
                if os.path.isfile(bin_path): break
            else:
                err_msg = ("Cannot find path of bin_name %s, neither in the build directory nor in PATH %s" %
                           (bin_name, paths))
                #warnings.warn(err_msg)
                bin_path = ""

        return bin_path

    def cygwin_path_of_bin(self, bin_name):
        """
        Mangle the name of the executable. Needed for Windows
        when we have to call an executable that is not located
        within the CYGWIN filesystem (aka $Win$ application).
        """
        path = self.path_of_bin(bin_name)
        if self.iscygwin(): path = self._cygwin_instdir + path
        return path

    def has_bin(self, bin_name, try_syspath=True):
        """True if binary bin_name is present in the build."""
        return os.path.isfile(self.path_of_bin(bin_name, try_syspath=try_syspath))

    def cygwin_path(self, path):
        apath = os.path.abspath(path)
        if self.iscygwin(): apath = self._cygwin_instdir + apath
        return apath

    def set_buildbot_builder(self, builder):
        """
        Set the name of the buildbot builder.
        Used to skip tests defining `exclude_builders` in the TEST_INFO_SECTION
        """
        self.buildbot_builder = builder


def parse_configh_file(fname):
    """
    Parse the configuration file config.h,
    Returns a list with the CCP variables that are #defined.

    Note:
      Not very robust. It does not handle instructions such as:

      #ifdef HAVE_FOO
      #  define HAVE_BAR 1
      #endif

    Handling this case would require a real preprocessing with CPP and then the parsing.
    Not easy to implement in a portable way especially on IBM machines with XLF.
    """
    with open(fname, "r") as fh:
        #defined_cppvars = []
        #for l in fh:
        #    l = l.lstrip()
        #    if l.startswith("#define "):
        #        tokens = l.split()
        #        varname = tokens[1]
        #        if varname.startswith("HAVE_") and len(tokens) >= 3:
        #            value = int(tokens[2])
        #            if value != 0: defined_cppvars.append(varname)

        defined_cppvars = {}
        for l in fh:
            l = l.lstrip()
            if l.startswith("#define "):
                tokens = l.split()
                varname = tokens[1]
                if len(tokens) >= 3:
                    value = tokens[2]
                    defined_cppvars[varname] = value

        return defined_cppvars


def input_file_has_vars(fname, ivars, comment="#", mode="any"):
    """
    Primitive parser that searches for the occurrence of input variables in the input file fname

    Args:
        fname:
            Input file
        ivars:
            dictionary whose keys are strings with the input variables to search.
            ivar[varname] can be either None or an integer
            if ivar[varname] is None, we have a match if varname is present
            if ivar[varname] is int, we have a match if varname is present and it has value int

        mode: "all" or "any"

        return:
            (bool, d)
            bool is True is the input file contains the specified variables
            d is a dictionary with the matching lines (empty dict if no occurence).
    """
    # This algorithm is not very robust as it assumes that the variable and the line
    # are placed on the same line.
    with open(fname, "r") as fh:
        lines = []
        for line in fh:
            line = line.lower().strip()
            idx = line.find(comment)
            if idx != -1: line = line[:idx]
            lines.append(line)

    matches = {}
    for k in ivars:
        matches[k] = []

    items = ivars.items()

    re_ivars = {}
    for varname in ivars:
        re_ivars[varname] = re.compile(varname + "\d*\s*(\d+)\s*")

    nfound = 0
    for line in lines:
        for varname, varvalue in items:
            re_match = re_ivars[varname].match(line)
            #print("match")
            if varvalue is None and varname in line:
                nfound += 1
                matches[varname].append(line)
            elif re_match:
                num = int(re_match.group(1))
                if num == int(varvalue):
                    #print line
                    matches[varname].append(line)
                    nfound += 1

    if nfound == 0:
        return False, {}

    if mode == "all":
        return all(bool(v) for v in matches.values()), matches
    elif mode == "any":
        return any(bool(v) for v in matches.values()), matches
    else:
        raise ValueError("Wrong mode %s" % mode)


class FldiffResult(object):
    """Store the results produced by fldiff.pl."""
    _attrbs = {
        "fname1": "first file provided to fldiff.",
        "fname2": "second file provided to fldiff.",
        "options": "options passed to fldiff.",
        "summary_line": "Summary given by fldiff.",
        "fatal_error": "True if file comparison cannot be done.",
        "ndiff_lines": "Number of different lines.",
        "abs_error": "Max absolute error.",
        "rel_error": "Max relative error.",
        "max_absdiff_ln": "Line number where the Max absolute error occurs.",
        "max_reldiff_ln": "Line number where the Max relative error occurs.",
    }

    def __init__(self, summary_line, err_msg, fname1, fname2, options):

        self.summary_line = summary_line.strip()
        self.err_msg = err_msg.strip()
        self.fname1 = fname1
        self.fname2 = fname2
        self.options = options

        self.fatal_error = False
        self.success = False

        if "fatal" in summary_line:
            self.fatal_error = True
        elif "no significant difference" in summary_line:
            self.success = True
            self.ndiff_lines = 0
            self.abs_error = 0.0
            self.rel_error = 0.0
        elif "different lines=" in summary_line:
            #Summary Case_84 : different lines= 5 , max abs_diff= 1.000e-03 (l.1003), max rel_diff= 3.704e-02 (l.1345)
            tokens = summary_line.split(",")
            for tok in tokens:
                if "different lines=" in tok:
                    self.ndiff_lines = int(tok.split("=")[1])
                if "max abs_diff=" in tok:
                    vals = tok.split("=")[1].split()
                    self.abs_error = float(vals[0])
                if "max rel_diff=" in tok:
                    vals = tok.split("=")[1].split()
                    self.rel_error = float(vals[0])
        else:
            err_msg = "Wrong summary_line: " + str(summary_line)
            #raise ValueError(err_msg)
            warnings.warn(err_msg)
            self.fatal_error = True

    @lazy__str__
    def __str__(self): pass

    def passed_within_tols(self, tolnlines, tolabs, tolrel):
        """
        Check if the test passed withing the specified tolerances.

        Returns:
            (isok, status, msg)
        """
        status = "succeeded"; msg = ""
        if self.fatal_error:
            status = "failed"
            msg = "fldiff.pl fatal error:\n" + self.err_msg
        elif self.success:
            msg = "succeeded"
        else:
            abs_error = self.abs_error
            rel_error = self.rel_error
            ndiff_lines = self.ndiff_lines
            status = "failed"; fact = 1.0

            locs = locals()
            if abs_error > tolabs * fact and rel_error < tolrel:
                msg = "failed: absolute error %(abs_error)s > %(tolabs)s" % locs
            elif rel_error > tolrel * fact and abs_error < tolabs:
                msg = "failed: relative error %(rel_error)s > %(tolrel)s" % locs
            elif ndiff_lines > tolnlines:
                msg = "failed: erroneous lines %(ndiff_lines)s > %(tolnlines)s" % locs
            elif abs_error > tolabs * fact and rel_error > tolrel * fact:
                msg = "failed: absolute error %(abs_error)s > %(tolabs)s, relative error %(rel_error)s > %(tolrel)s" % locs
            # FIXME passed or failed?
            elif abs_error > tolabs:
                msg = "within 1.5 of tolerance (absolute error %(abs_error)s, accepted %(tolabs)s )" % locs
            elif rel_error > tolrel:
                msg = "within 1.5 of tolerance (relative error %(rel_error)s, accepted %(tolrel)s )" % locs
            else:
                status = "passed"
                msg = "passed: absolute error %(abs_error)s < %(tolabs)s, relative error %(rel_error)s < %(tolrel)s" % locs

        isok = status in ["passed", "succeeded"]

        #if not self.success:
        # Add the name of the file.
        msg += " [file=%s]" % os.path.basename(self.fname1)

        return isok, status, msg


def wrap_fldiff(fldiff_path, fname1, fname2, opts=None, label=None, timebomb=None, out_filobj=sys.stdout):
    """
    Wraps fldiff.pl script, returns (fld_result, got_summary)

    fld_result is a FldiffResult instance, got_summary is set to False if fldiff.pl didn't return any final summary

    Usage: fldiff [-context] [ -ignore | -include ] [ -ignoreP | -includeP ] [ -easy | -medium | -ridiculous ] file1 file2 [label]
    """
    # Default options for fldiff script.
    fld_options = "-ignore -ignoreP"
    if opts: fld_options = " ".join([fld_options] + [o for o in opts])
    fld_options = [s for s in fld_options.split()]

    if label is None: label = ""

    args = ["perl", fldiff_path] + fld_options + [fname1, fname2, label]
    cmd_str = " ".join(args)

    logger.info("about to execute %s" % cmd_str)

    if True or timebomb is None:
        if py2:
            p = Popen(cmd_str, shell=True, stdout=PIPE, stderr=PIPE)
        else:
            p = Popen(cmd_str, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        stdout_data, stderr_data = p.communicate()
        ret_code = p.returncode
        #ret_code = p.wait()
    else:
        p, ret_code = timebomb.run(cmd_str, shell=True, stdout=PIPE, stderr=PIPE)

    # fldiff returns this value when some difference is found.
    # perl programmers have a different understanding of exit_status!
    MAGIC_FLDEXIT = 4

    err_msg = ""
    if ret_code not in [0, MAGIC_FLDEXIT]:
        #err_msg = p.stderr.read()
        err_msg = stderr_data

    lines = stdout_data.splitlines(True)
    #lines = [str(s) for s in stdout_data.splitlines(True)]
    #lines = p.stdout.readlines()

    if out_filobj and not hasattr(out_filobj, "writelines"):
        # Assume string
        lazy_writelines(out_filobj, lines)
    else:
        #print(out_filobj)
        out_filobj.writelines(lines)

    # Parse the last line.
    # NOTE:
    # on woopy fldiff returns to the parent process without producing
    # any output. In this case, we set got_summary to False so that
    # the caller can make another attempt.
    got_summary = True

    try:
        summary_line = lines[-1]
    except IndexError:
        got_summary = False
        try:
            logger.critical("Trying to kill fldiff process, cmd %s" % cmd_str)
            p.kill()
        except Exception as exc:
            logger.critical("p.kill failed with exc %s" % str(exc))
            pass
        summary_line = "fatal error: no summary line received from fldiff"

    return FldiffResult(summary_line, err_msg, fname1, fname2, fld_options), got_summary


def make_abitest_from_input(inp_fname, abenv, keywords=None, need_cpp_vars=None, with_np=1):
    """
    Factory function to generate a Test object from the input file inp_fname
    """
    inp_fname = os.path.abspath(inp_fname)

    try: # Assumes some_path/Input/t30.in
        inpdir_path, x = os.path.split(inp_fname)
    except:
        raise ValueError("%s is not a valid path" % inp_fname)

    parser = AbinitTestInfoParser(inp_fname)

    nprocs_to_test = parser.nprocs_to_test
    ntests = len(nprocs_to_test)

    if ntests == 0:
        nprocs_to_test = [1]
        ntests = 1

    test_info = parser.generate_testinfo_nprocs(with_np)

    # Add global cpp variables.
    test_info.add_cpp_vars(need_cpp_vars)

    # Add global keywords.
    test_info.add_keywords(keywords)

    # Single test with np processors.
    # Istanciate the appropriate subclass depending on the name of the executable. Default is BaseTest.
    cls = exec2class(test_info.executable)

    return cls(test_info, abenv)


def make_abitests_from_inputs(input_fnames, abenv, keywords=None, need_cpp_vars=None):
    """
    Factory function. Return a list of tests generated from the TEST_INFO section reported
    in the input files inp_fnames.
    """
    if is_string(input_fnames):
        input_fnames = [input_fnames]

    inp_fnames = [os.path.abspath(p) for p in input_fnames]

    out_tests = []

    while True:
        try:
            inp_fname = inp_fnames.pop(0)
        except IndexError:
            break

        try:
            # Assumes some_path/Input/t30.in
            inpdir_path, x = os.path.split(inp_fname)
        except:
            raise ValueError("%s is not a valid path" % inp_fname)

        #print("inp_fname", inp_fname)
        parser = AbinitTestInfoParser(inp_fname)
        nprocs_to_test = parser.nprocs_to_test

        if len(nprocs_to_test) == 0:
            nprocs_to_test = [1]

        if not parser.is_testchain:
            # No dependency --> generate a list of test by changing the number np of MPI processors.
            for np in nprocs_to_test:
                test_info = parser.generate_testinfo_nprocs(np)

                test_info.add_cpp_vars(need_cpp_vars) # Add global cpp variables.
                test_info.add_keywords(keywords)      # Add global keywords.

                # Istanciate the appropriate subclass depending on the name of the executable. Default is BaseTest.
                cls = exec2class(test_info.executable)
                out_tests.append(cls(test_info, abenv))

        else:
            logger.info("got chain input %s" % inp_fname)
            #print(parser.chain_inputs())

            # Build the test chain with np nprocessors.
            for np in nprocs_to_test:
                tchain_list = []
                for cht_fname in parser.chain_inputs():
                    t = make_abitest_from_input(cht_fname, abenv, keywords=keywords, need_cpp_vars=need_cpp_vars, with_np=np)
                    tchain_list.append(t)

                if not tchain_list:
                    raise RuntimeError("tchain_list is empty, inp_fname %s" % inp_fname)

                out_tests.append(ChainOfTests(tchain_list))

            # Remove the input files of the chain
            for s in parser.chain_inputs()[1:]:
                try:
                    idx = inp_fnames.index(s)
                except ValueError:
                    raise RuntimeError("%s not found in inp_fnames" % inp_fnames)

                inp_fnames.pop(idx)

    return out_tests


class Status(int):
    """
    This object is an integer representing the status of the `Test`.

    Statuses are ordered, negative values are used for positive outcomes,
    positive values for failures.
    """
    # Possible status of the node.
    _STATUS2STR = OrderedDict([
        (-3, "Skipped"),         # Test has been skipped because test requirements are not fulfilled
        (-2, "Succeeded"),       # fldiff returned succeeded
        (-1, "Passed"),          # fldiff returned passed
        (0, "None"),             # Initial status of the test.
        (1, "FileDifferError"),  # File comparison could not be performed but the calculation terminated
                                 # (e.g. different number of lines in ref and out files)
        (2, "NumericalError"),   # File comparison detected too large numerical errors.
        (3, "ExecutionError"),   # Run didn't complete due to some error in the code e.g. segmentation fault
        (4, "PythonError"),      # A python exception was raised in the driver code.
    ])

    def __repr__(self):
        return "<%s: %s, at %s>" % (self.__class__.__name__, str(self), id(self))

    def __str__(self):
        """String representation."""
        return self._STATUS2STR[self]

    @classmethod
    def from_string(cls, s):
        """Return a `Status` instance from its string representation."""
        for num, text in cls._STATUS2STR.items():
            if text == s:
                return cls(num)
        else:
            raise ValueError("Wrong string %s" % s)

    @property
    def is_problematic(self):
        """True if test was not successful."""
        return self > 0

    @property
    def info(self):
        """Human-readable string with info on the outcome of the test."""
        try:
            return self._info
        except AttributeError:
            return "None"

    def set_info(self, info):
        """info setter."""
        self._info = info


class BaseTestError(Exception):
    """Base Error class raised by Test objects"""


class BaseTest(object):
    """
    Base class describing a single test. Tests associated to other executables should
    sublcass BaseTest and redefine the method make_stdin.
    Then change exec2cls so that the appropriate instance is returned.
    """
    Error = BaseTestError

    # Possible status of the test.
    _possible_status = ["failed", "passed", "succeeded", "skipped", "disabled"]

    #S_SKIPPED = Status(-3)
    #S_SUCCEDED = Status(-2)
    #S_PASSED = Status(-1)
    #S_NODE = Status(0)
    #S_DIFF_ERROR = Status(2)
    #S_NUM_ERROR = Status(2)
    #S_EXEC_ERROR = Status(3)
    #S_PY_ERROR = Status(4)

    #ALL_STATUS = [
    #    S_SKIPPED,
    #    S_SUCCEDED,
    #    S_PASSED,
    #    S_NODE,
    #    S_NUM_ERROR,
    #    S_EXEC_ERROR,
    #    S_PY_ERROR,
    #]

    def __init__(self, test_info, abenv):
        logger.info("Initializing BaseTest from inp_fname: ", test_info.inp_fname)

        self.inp_fname = os.path.abspath(test_info.inp_fname)
        self.abenv = abenv
        self.id = test_info.make_test_id() # The test identifier (takes into account the multi_parallel case)
        self.nprocs = 1  # Start with 1 MPI process.

        # FIXME Assumes inp_fname is in the form tests/suite_name/Input/name.in
        suite_name = os.path.dirname(self.inp_fname)
        suite_name = os.path.dirname(suite_name)

        self.suite_name = os.path.basename(suite_name)
        self.ref_dir = abenv.apath_of("tests", suite_name, "Refs")
        self.inp_dir = abenv.apath_of("tests", suite_name, "Input")

        self._executed = False
        self._status = None
        if os.path.basename(self.inp_fname).startswith("-"):
            self._status = "disabled"

        # Initial list of local files that should not be removed.
        self._files_to_keep = []

        # Default values.
        self.make_html_diff = 0   # 0 => Do not produce diff files in HTML format
                                  # 1 => Produced HTML diff but only if test failed
                                  # 2 => Produce HTML diff independently of the final status

        self.sub_timeout = 30     # Timeout for subprocesses (in seconds)

        self.erase_files = 2      # 0 => Keep all files.
                                  # 1 => Remove files but only if the test passes or succeeds
                                  # 2 => Remove files even when the test fail.

        # Incorporate the attributes of test_info in self.
        err_msg = ""
        for k in test_info.__dict__:
            if k in self.__dict__ and test_info.__dict__[k] != self.__dict__[k]:
                err_msg += "Cannot overwrite key %s\n" % k
                #print(test_info.__dict__[k],  self.__dict__[k])

        if err_msg:
            raise self.Error(err_msg)

        self.__dict__.update(test_info.__dict__)

        # Save authors' second names to speed up the search.
        # Well, let's hope that we don't have authors with the same second name!
        second_names = []
        for string in self.authors:
            idx = string.rfind(".")
            f, s = ("", string)
            if idx != -1:
                try:
                    f, s = string[:idx+1], string[idx+2:]
                except IndexError:
                    raise ValueError("Wrong author(s) name")

            if not f and s and s != "Unknown":
                print("author(s) first name is missing in file %s, string = %s " % (self.full_id, string))

            second_names.append(s)

        self._authors_snames = set(second_names)

    def __repr__(self):
        return self.full_id

    def __str__(self):
        return repr(self)

    #@lazy__str__
    #def __str__(self): pass

    def stdin_readlines(self):
        return lazy_readlines(self.stdin_fname)

    def stdin_read(self):
        return lazy_read(self.stdin_fname)

    def stdout_readlines(self):
        return lazy_readlines(self.stdout_fname)

    def stdout_read(self):
        return lazy_read(self.stdout_fname)

    def stderr_readlines(self):
        return lazy_readlines(self.stderr_fname)

    def stderr_read(self):
        return lazy_read(self.stderr_fname)

    @property
    def has_empty_stderr(self):
        return not bool(self.stderr_read())

    @property
    def full_id(self):
        """Full identifier of the test."""
        return "[%s][%s][np=%s]" % (self.suite_name, self.id, self.nprocs)

    @property
    def bin_path(self):
        """The absolute path of the executable needed to run the test."""
        return self.build_env.path_of_bin(self.executable)

    @property
    def cygwin_bin_path(self):
        return self.build_env.cygwin_path_of_bin(self.executable)

    def cygwin_path(self, path):
        return self.build_env.cygwin_path(path)

    def cpkl_dump(self, protocol=-1):
        """Save the instance in a pickle file"""
        self.cpkl_fname = os.path.join(self.workdir, self.id + ".cpkl")

        with open(self.cpkl_fname, "wb") as fh:
            pickle.dump(self, fh, protocol=protocol)
            self.keep_files(self.cpkl_fname)

    def has_keywords(self, keywords, mode="any"):
        """
        True if test has keywords
        mode == "all" --> check if all keywords are present
        mode == "any" --> check if at least one keyword is present
        """
        if mode == "all":
            return set(keywords).issubset(self.keywords)
        elif mode == "any":
            return set(keywords).intersection(self.keywords)
        else:
            raise ValueError("wrong mode %s" % mode)

    def has_authors(self, authors, mode="any"):
        """
        True if test has authors

        mode == "all" --> check if all authors are present
        mode == "any" --> check if at least one author is present
        """
        if mode == "all":
            return set(authors).issubset(self._authors_snames)
        elif mode == "any":
            return set(authors).intersection(self._authors_snames)
        else:
            raise ValueError("wrong mode %s" % mode)

    def get_varname_set(self):
        """
        Return set of variables used by this test.
        Mainly used to check if all variables in the doc are documented/tested.

        .. note:

            Dataset index (if any) is removed.
        """
        # See abio.abivars.AbinitInputParser
        import io
        lines = []
        with io.open(self.inp_fname, "rt", encoding="utf-8") as fh:
            for line in fh:
                line.strip()
                # Remove comments from lines.
                i = line.find("#")
                if i != -1: line = line[:i]
                i = line.find("!")
                if i != -1: line = line[:i]
                if line: lines.append(line)

        vnames = []
        # 1) Build string of the form "var1 value1 var2 value2"
        tokens = " ".join(lines).split()
        for pos, tok in enumerate(tokens):
            if tok[0].isalpha():
                # TODO
                # Either new variable, string defining the unit or operator e.g. sqrt
                #if is_abiunit(tok) or tok in ABI_OPERATORS or "?" in tok:
                #    continue

                # Have new variable
                if tok[-1].isdigit(): # and "?" not in tok:
                    # Handle dataset index.
                    l = []
                    for i, c in enumerate(tok[::-1]):
                        if c.isalpha(): break
                        #l.append(c)
                    else:
                        raise ValueError("Cannot find dataset index in token: %s" % tok)
                    tok = tok[:len(tok) - i]
                    #l.reverse()
                    #print("tok", tok, l)
                    #tok = l
                vnames.append(tok)

        #print(vnames)
        return set(v.lower() for v in vnames)

    def has_variables(self, ivars, mode="any"):
        """True if test has the input variables ivars (dict {varname:varvalue})"""
        found, d = input_file_has_vars(self.inp_fname, ivars, mode=mode)
        return found

    def edit_input(self, editor=None):
        """
        Call editor to edit the input file of the test.
        A default editor is provided if editor is None (use $EDITOR shell variable)
        """
        if editor is None: editor = Editor()
        try:
            editor.edit_file(self.inp_fname)
        except:
            raise

    def listoftests(self, width=100, html=True, abslink=True):
        string = self.description.lstrip()
        if self.references:
            string += "References:\n" + "\n".join(self.references)
        string = textwrap.dedent(string)
        string = textwrap.fill(string, width=width)
        if not html:
            return self.full_id + ":\n" + string
        else:
            if abslink:
                link = html_link(self.full_id, self.inp_fname)
            else:
                # Use relative path so that we can upload the HTML file on
                # the buildbot master and browse the pages.
                link = html_link(self.full_id, os.path.basename(self.inp_fname))
            string = link + "<br>" + string.replace("\n","<br>") + "\n"
        return string

    def make_stdin(self):
        """
        Generate the standard input of the test.
        The base implementation writes the content of inp_fname to stdin.
        Subclasses should redefine this method according to their needs.
        """
        t_stdin = StringIO()
        with open(self.inp_fname, "r") as fh:
            t_stdin.writelines(fh)

        return t_stdin.getvalue()

    def get_extra_inputs(self):
        """Copy extra inputs from inp_dir to workdir."""
        # First copy the main input file (useful for debugging the test)
        # Avoid raising exceptions as python threads do not handle them correctly.
        try:
            src = self.inp_fname
            dest = os.path.join(self.workdir, os.path.basename(self.inp_fname))
            shutil.copy(src, dest)
            self.keep_files(dest)  # Do not remove it after the test.
        except:
            self.exceptions.append(self.Error("copying %s => %s" % (src,dest)))

        for extra in self.extra_inputs:
            src = os.path.join(self.inp_dir, extra)
            dest = os.path.join(self.workdir, extra)

            if not os.path.isfile(src):
                self.exceptions.append(self.Error("%s: no such file" % src) )
                continue

            shutil.copy(src, dest)
            if dest.endswith(".gz"):  # Decompress the file
                unzip(dest)
                dest = dest[:-3]
            #self.keep_files(dest)  # Do not remove dest after the test.

    @property
    def inputs_used(self):
        """List with the input files used by the test."""
        inputs = [self.inp_fname] + [os.path.join(self.inp_dir, f) for f in self.extra_inputs]
        #
        # Add files appearing in the shell sections.
        for cmd_str in (self.pre_commands + self.post_commands):
            if cmd_str.startswith("iw_"):
                tokens = cmd_str.split()
                inp = os.path.join(self.inp_dir, tokens[1])
                inputs.append(inp)

        return inputs

    @property
    def status(self):
        """The status of the test"""
        if self._status in ["disabled", "skipped", "failed"]: return self._status
        all_fldstats = [f.fld_status for f in self.files_to_test]
        if "failed" in all_fldstats: return "failed"
        if "passed" in all_fldstats: return "passed"
        assert all([s == "succeeded" for s in all_fldstats])

        return "succeeded"

    @property
    def isok(self):
        """Return true if test is OK (test passed and not python exceptions."""
        return self.fld_isok and not self.exceptions

    @property
    def files_to_keep(self):
        """List with the files that should not be erased once the test completed"""
        return self._files_to_keep

    def keep_files(self, files):
        """Add files to the list of paths that should not be erased"""
        if is_string(files):
            self._files_to_keep.append(files)
        else:
            self._files_to_keep.extend(files)

    def compute_nprocs(self, build_env, nprocs, runmode):
        """
        Compute the number of MPI processes that can be used for the test from the initial guess nprocs

        Return: (nprocs, string)

        where nprocs = 0 if the test cannot be executed.
        string contains a human-readable message explaining the reason why the test will be skipped.

        A test cannot be executed if:

          1) It requires CPP variables that are not defined in the build.
          2) The user asks for more MPI nodes than max_nprocs (this value is reported in the TEST_INFO section).
          3) We have a multiparallel test (e.g. paral/tA.in) and nprocs is not in in nprocs_to_test
          4) nprocs is in exclude_nprocs
        """
        # !HAVE_FOO --> HAVE_FOO should not be present.
        errors = []
        eapp = errors.append
        for var in self.need_cpp_vars:
            if not var.startswith("!") and var not in build_env.defined_cppvars:
                eapp("Build environment does not define the CPP variable %s" % var)
            elif var[1:] in build_env.defined_cppvars:
                eapp("Build environment defines the CPP variable %s" % var[1:])

        # Remove this check to run the entire test suite in parallel
        #runmode ="dynamic"

        if runmode == "static":
            if nprocs > self.max_nprocs:
                eapp("nprocs: %s > max_nprocs: %s" % (nprocs, self.max_nprocs))

        elif runmode == "dynamic":
            # Will select the minimum between max_nprocs and nprocs
            pass

        else:
            raise ValueError("Wrong runmode %s" % runmode)

        if self.nprocs_to_test and nprocs != self.nprocs_to_test[0]:
            eapp("nprocs: %s != nprocs_to_test: %s" % (nprocs, self.nprocs_to_test[0]))

        if nprocs in self.exclude_nprocs:
            eapp("nprocs: %s in exclude_nprocs: %s" % (nprocs, self.exclude_nprocs))

        err_msg = "\n".join(errors)
        if err_msg:
            real_nprocs = 0
        else:
            real_nprocs = min(self.max_nprocs, nprocs)

        #if err_msg: print(err_msg)
        return real_nprocs, err_msg

    def skip_host(self):
        """
        Return True if the test should be skipped since we are running on a banned host.
        """
        compilers, slaves = [], []

        for s in self.exclude_hosts:
            compiler, host = None, s
            if "@" in s:
                compiler, host = s.split("@")
            else:
                # TODO: validate TEST_INFO at the level of the parser.
                warnings.warn("Wrong string %s in exclude_hosts" % s)

            compilers.append(compiler)
            slaves.append(host)

        # Find the slave and compare the name of the compiler.
        try:
            idx = slaves.index(self.build_env.hostname)
        except ValueError:
            return False

        return compilers[idx] == self.build_env.fortran_compiler.name

    def skip_buildbot_builder(self):
        """
        Return True if the test should be skipped since we are running on a banned builder.
        """
        if not hasattr(self.build_env, "buildbot_builder"): return False
        for builder in self.exclude_builders:
            if builder == self.build_env.buildbot_builder: return True
        return False

    def run(self, build_env, runner, workdir, nprocs=1, runmode="static", **kwargs):
        """
        Run the test with nprocs MPI nodes in the build environment build_env using the `JobRunner` runner.
        Results are produced in directory workdir. kwargs is used to pass additional options

        ================  ====================================================================
        kwargs            Meaning
        ================  ====================================================================
        pedantic           Mark tests as failed if stderr is not empty.
        erase_files        0 => Keep all files produced by the test
                           1 => Remove files but only if the test passed or succeeded.
                           2 => Remove files even if the test failed.
                           default=2
        make_html_diff     True to produce diff in HTML format. Default: False.
        sub_timeout        Timeout for subprocesses.
        abimem_check       True if abimem.mocc files should be analyzes for possible errors.
                           Requires HAVE_MEM_PROFILE and `call abimem_init(2)` in main.
                           Default: False
        etsf_check         True if netcdf files should be validated. Requires netcdf4.
                           Default: False
        ================  ====================================================================

        .. warning:
            This method must be thread-safe, DO NOT change build_env or runner.
        """
        import copy
        runner = copy.deepcopy(runner)
        start_time = time.time()

        workdir = os.path.abspath(workdir)
        if not os.path.exists(workdir): os.mkdir(workdir)
        self.workdir = workdir

        self.build_env = build_env

        self.exceptions = []
        self.fld_isok = True  # False if at least one file comparison fails.

        # Extract options from kwargs
        self.pedantic = kwargs.get("pedantic", False)
        self.erase_files = kwargs.get("erase_files", self.erase_files)
        self.make_html_diff = kwargs.get("make_html_diff", self.make_html_diff)
        self.sub_timeout = kwargs.get("sub_timeout", self.sub_timeout)


        timeout = self.sub_timeout
        if self.build_env.has_bin("timeout") and timeout > 0.0:
            exec_path = self.build_env.path_of_bin("timeout")
            self.timebomb = TimeBomb(timeout, delay=0.05, exec_path = exec_path)
        else:
            self.timebomb = TimeBomb(timeout, delay=0.05)

        str_colorizer = StringColorizer(sys.stdout)

        #status2txtcolor = {
        #    "succeeded": lambda string: str_colorizer(string, "green"),
        #    "passed": lambda string: str_colorizer(string, "blue"),
        #    "failed": lambda string: str_colorizer(string, "red"),
        #    "disabled": lambda string: str_colorizer(string, "cyan"),
        #    "skipped": lambda string: str_colorizer(string, "cyan"),
        #}

        status2txtcolor = {
            "succeeded": "green",
            "passed": "blue",
            "failed": "red",
            "disabled": "cyan",
            "skipped": "cyan",
        }

        # Check whether the test can be executed.
        can_run = True
        if self._status == "disabled":
            msg = self.full_id + ": Disabled"
            can_run = False
            cprint(msg, status2txtcolor[self._status])

        # Here we get the number of MPI nodes for test.
        self.nprocs, self.skip_msg = self.compute_nprocs(self.build_env, nprocs, runmode=runmode)

        if self.skip_msg:
            self._status = "skipped"
            msg = self.full_id + ": Skipped."
            cprint(msg, status2txtcolor[self._status])
            for l in self.skip_msg.splitlines():
                cprint("\t" + l, status2txtcolor[self._status])
            print()
            can_run = False

        if self.skip_host():
            self._status = "skipped"
            msg = self.full_id + ": Skipped: this hostname has been excluded."
            cprint(msg, status2txtcolor[self._status])
            can_run = False

        if self.skip_buildbot_builder():
            self._status = "skipped"
            msg = self.full_id + ": Skipped: this buildbot builder has been excluded."
            cprint(msg, status2txtcolor[self._status])
            can_run = False

        self.run_etime = 0.0

        if can_run:
            # Execute pre_commands in workdir.
            rshell = RestrictedShell(self.inp_dir, self.workdir, self.abenv.psps_dir)

            for cmd_str in self.pre_commands:
                rshell.execute(cmd_str)

            if rshell.exceptions:
                self.exceptions.extend(rshell.exceptions)
                rshell.empty_exceptions()

            # Copy extra inputs in workdir (if any).
            self.get_extra_inputs()

            # Create stdin file in the workdir.
            self.stdin_fname = os.path.join(workdir, self.id + ".stdin")
            self.stdout_fname = os.path.join(workdir, self.id + ".stdout")
            self.stderr_fname = os.path.join(workdir, self.id + ".stderr")

            self.keep_files([self.stdin_fname, self.stdout_fname, self.stderr_fname])

            # Create input file.
            t_stdin = self.make_stdin()
            with open(self.stdin_fname, "w") as fh:
                fh.writelines(t_stdin)

            # Run the code (run_etime is the wall time spent to execute the test)
            if runner.has_mpirun:
                bin_path = self.cygwin_bin_path
            else:
                bin_path = self.bin_path

            self.run_etime = runner.run(self.nprocs, bin_path,
                                        self.stdin_fname, self.stdout_fname, self.stderr_fname,
                                        cwd=workdir)

            # Save exceptions (if any).
            if runner.exceptions:
                self.exceptions.extend(runner.exceptions)
                if not self.expected_failure:
                    for exc in runner.exceptions: print(exc)

            # Execute post_commands in workdir.
            for cmd_str in self.post_commands:
                rshell.execute(cmd_str)

            # Save exceptions (if any).
            if rshell.exceptions:
                self.exceptions.extend(rshell.exceptions)
                rshell.empty_exceptions()

            # Check final results:
            # 1) use fldiff to compare ref and output files.
            # 2) fldiff stdout is redirected to fldiff_fname.
            for f in self.files_to_test:
                fldiff_fname = os.path.join(self.workdir, f.name + ".fldiff")
                self.keep_files(fldiff_fname)

                with open(fldiff_fname,"w") as fh:
                    f.fldiff_fname = fldiff_fname

                    isok, status, msg = f.compare(self.abenv.fldiff_path, self.ref_dir, self.workdir,
                                                  timebomb=self.timebomb, outf=fh)

                self.keep_files(os.path.join(self.workdir, f.name))
                self.fld_isok = self.fld_isok and isok

                msg = ": ".join([self.full_id, msg])
                cprint(msg, status2txtcolor[status])

            # Check if the test is expected to fail.
            if runner.retcode != 0 and not self.expected_failure:
                self._status = "failed"
                msg = (self.full_id + "Test was not expected to fail but subprocesses returned %s" % runner.retcode)
                cprint(msg, status2txtcolor["failed"])

            # If pedantic, stderr must be empty unless the test is expected to fail!
            if self.pedantic and not self.expected_failure:
                try:
                    errout = self.stderr_read()
                    if errout:
                        # TODO: Not very clean, I should introduce a new status and a setter method.
                        self._status = "failed"
                except Exception as exc:
                    self.exceptions.append(exc)

            # Check stderr for presence of valgrind errors.
            if runner.has_valgrind:
                try:
                    # Build a parser from the command line options and parse the stderr.
                    parser = runner.build_valgrind_parser()
                    parser.parse(self.stderr_fname)

                    if parser.error_report:
                        # TODO: Not very clean, I should introduce a new status and a setter method.
                        self._status = "failed"
                        msg = " ".join([self.full_id, "VALGRIND ERROR:", parser.error_report])
                        cprint(msg, status2txtcolor["failed"])

                except Exception as exc:
                    # Py threads do not like exceptions.
                    # Store the exception and continue.
                    self.exceptions.append(exc)

            if self.status == "failed":
                # Print the first line of the stderr if it's not empty.
                # Look also for the MPIABORTFILE
                try:
                    errout= self.stderr_read()
                    if errout:
                        cprint(errout, status2txtcolor["failed"])

                    # Extract YAML error message from ABORTFILE or stdout.
                    abort_file = os.path.join(self.workdir, "__ABI_MPIABORTFILE__")
                    if os.path.exists(abort_file):
                        f = open(abort_file, "r")
                        print( 12*"=" + " ABI_MPIABORTFILE " + 12*"=")
                        cprint(f.read(), status2txtcolor["failed"])
                        f.close()
                    else:
                        yamlerr = read_yaml_errmsg(self.stdout_fname)
                        if yamlerr:
                            print("YAML Error found in the stdout of", self)
                            cprint(yamlerr, status2txtcolor["failed"])
                        else:
                            print("No YAML Error found in", self)

                except Exception as exc:
                    self.exceptions.append(exc)

            if kwargs.get("abimem_check", False):
                paths = [os.path.join(self.workdir, f) for f in os.listdir(self.workdir)
                         if f.startswith("abimem") and f.endswith(".mocc")]
                print("Found %s abimem files" % len(paths))
                #abimem_retcode = 0
                for path in paths:
                    parser = AbimemParser(path)
                    parser.find_memleaks()
                    #if rc: parser.show_errors()
                    #abimem_retcode += rc

            #if False and kwargs.get("etsf_check", False):
            if kwargs.get("etsf_check", False):
		# Mark the test as failed and create a custom Exception
		# developers will have to inspect the xreport file for the full list of errors.
                try:
                    from . import etsf_specs as etsf
                except ImportError:
                    etsf = None

                errmsg = ""
                if etsf is None:
                    errmsg ="etsf_check is activated but netcdf4 module is not available"
                    nc_retcode = 1

                else:
                    nc_retcode = 0
                    all_errors = []
                    for p in os.listdir(self.workdir):
                        if not p.endswith(".nc"): continue
                        path = os.path.join(self.workdir, p)
                        elist = []
                        #elist += etsf.validate_vars(path))
                        elist += etsf.validate_ncfile(path)

                        if elist:
                            all_errors.append(elist)
                            cprint("%s [FAILED]" % p, "red")
                        else:
                            cprint("%s [OK]" % p, "green")

                    nc_retcode = len(all_errors)

                    if nc_retcode != 0:
                        errmsg = ("Setting status to failed because nc_retcode=%s\n"
                                  "The netcdf files produced by this tests either is not consistent with the etsf specs.\n"
                                  "or it has not been registered in ~abinit/tests/pymods/etsf_specs.py\n"
                                  "Please, control the errors messages in the xreport file produced by buildbot."
                                   % nc_retcode)

                if nc_retcode != 0:
                    # TODO: Not very clean, I should introduce a new status and a setter method.
                    self._status = "failed"
                    # Store the exception and continue.
                    self.exceptions.append(Exception(errmsg))
                    print(errmsg)
                else:
                    cprint("netcdf validation [OK]", "green")

        self._executed = True
        self.tot_etime = time.time() - start_time

    @property
    def executed(self):
        return self._executed

    def clean_workdir(self, other_test_files=None):
        """Remove the files produced in self.workdir."""
        assert self._executed
        if not os.path.exists(self.workdir) or self.erase_files == 0: return

        save_files = self._files_to_keep[:]
        if other_test_files is not None: save_files += other_test_files

        # Add harcoded list of files
        hard_files = ["perf.data", "__ABI_MPIABORTFILE__"]
        save_files += [os.path.join(self.workdir, f) for f in hard_files]

        # List of file extensions to be preserved.
        keep_exts = [".flun", ".mocc"]

        if (self.erase_files == 1 and self.isok) or self.erase_files == 2:
            entries = [os.path.join(self.workdir, e) for e in os.listdir(self.workdir)]
            for entry in entries:
                if entry in save_files: continue
                _, ext = os.path.splitext(entry)
                if ext in keep_exts: continue
                if os.path.isfile(entry):
                    try:
                        os.remove(entry)
                    except OSError:
                        pass
                else:
                    raise NotImplementedError("Found directory: %s in workdir!!" % entry)

    def patch(self, patcher=None):
        """
        Patch the output files of the test with the specified patcher.
        A default patcher is provided if patcher is None (use $PATCHER shell variable)
        """
        assert self._executed
        for f in self.files_to_test:
            ref_fname = os.path.abspath(os.path.join(self.ref_dir, f.name))
            out_fname = os.path.abspath(os.path.join(self.workdir,f.name) )
            raise NotImplementedError("patcher should be tested")
            from tests.pymods import Patcher
            Patcher(patcher).patch(out_fname, ref_fname)

    def make_html_diff_files(self):
        """Generate and write diff files in HTML format."""
        assert self._executed
        if (self.make_html_diff == 0 or
            self._status in ["disabled", "skipped"]): return

        diffpy = self.abenv.apath_of("tests", "pymods", "diff.py")

        for f in self.files_to_test:
            if f.fld_isok and self.make_html_diff == 1:
                continue

            ref_fname = os.path.abspath(os.path.join(self.ref_dir, f.name))

            if not os.path.isfile(ref_fname) and ref_fname.endswith(".stdout"):
                ref_fname = ref_fname[:-7] + ".out"  # FIXME Hack due to the stdout-out ambiguity

            out_fname = os.path.abspath(os.path.join(self.workdir,f.name))

            # Check whether output and ref file exist.
            out_exists = os.path.isfile(out_fname)
            ref_exists = os.path.isfile(ref_fname)

            hdiff_fname = os.path.abspath(os.path.join(self.workdir, f.name + ".diff.html"))

            f.hdiff_fname = hdiff_fname

            x, ext = os.path.splitext(f.name)
            safe_hdiff = ext in [".out", ".stdout"] # Create HTML diff file only for these files

            if ref_exists and out_exists and safe_hdiff:
                out_opt = "-u"
                #out_opt = "-t"   # For simple HTML table. (can get stuck)
                #args = ["python", diffpy, out_opt, "-f " + hdiff_fname, out_fname, ref_fname ]
                args = [diffpy, out_opt, "-f " + hdiff_fname, out_fname, ref_fname ]
                cmd = " ".join(args)
                #print("Diff", cmd)

                p, ret_code = self.timebomb.run(cmd, shell=True, cwd=self.workdir)

                if ret_code != 0:
                    err_msg = "Timeout error (%s s) while executing %s, retcode = %s" % (
                      self.timebomb.timeout, str(args), ret_code)
                    self.exceptions.append(self.Error(err_msg))
                else:
                    self.keep_files(hdiff_fname)

    def make_txt_diff_files(self):
        """Generate and write diff files in txt format."""
        assert self._executed
        if self._status in ["disabled", "skipped"]:
            return

        #print(self._status)
        #if self._status not in ["failed", "passed"]:
        #  return
        diffpy = self.abenv.apath_of("tests", "pymods", "diff.py")

        for f in self.files_to_test:
            #print(f, f.fld_isok)
            if f.fld_isok:
                continue

            ref_fname = os.path.abspath(os.path.join(self.ref_dir, f.name))

            if not os.path.isfile(ref_fname) and ref_fname.endswith(".stdout"):
                ref_fname = ref_fname[:-7] + ".out"  # FIXME Hack due to the stdout-out ambiguity

            out_fname = os.path.abspath(os.path.join(self.workdir, f.name))

            # Check whether output and ref file exist.
            out_exists = os.path.isfile(out_fname)
            ref_exists = os.path.isfile(ref_fname)

            diff_fname = os.path.abspath(os.path.join(self.workdir, f.name + ".diff"))

            f.diff_fname = diff_fname

            x, ext = os.path.splitext(f.name)

            if ref_exists and out_exists:
                # n is for ndiff format, c for context, u for unified
                #for out_opt in ["-n", "-c"]:
                #out_opt = "-n"
                #out_opt = "-c"
                out_opt = "-u"
                args = [diffpy, out_opt, "-f " + diff_fname, out_fname, ref_fname]
                cmd = " ".join(args)

                (p, ret_code) = self.timebomb.run(cmd, shell=True, cwd=self.workdir)

                if ret_code != 0:
                    err_msg = "Timeout error (%s s) while executing %s, retcode = %s" % (
                      self.timebomb.timeout, str(args), ret_code)
                    self.exceptions.append(self.Error(err_msg))
                else:
                    self.keep_files(diff_fname)

    def write_html_report(self, fh=None, oc="oc"):
        """Write the HTML file summarizing the results of the test."""
        assert self._executed

        close_fh = False
        if fh is None:
            close_fh = True
            html_report = os.path.join(self.workdir, "test_report.html")
            fh = open(html_report, "w")

        self.keep_files(fh.name)

        self.make_html_diff_files()
        self.make_txt_diff_files()

        # Try to read stdout, stderr and the abort_file produced by Abinit in parallel
        # Ignore errors (fock takes years to flush the stdout)
        #stdout_text, stderr_text = 2*("",)
        nlast = 120
        stderr_text, stdout_text, abiabort_text = 3 * (" ",)
        abort_file = find_abortfile(self.workdir)
        #self.fld_isok = False
        errinfo_text = " "
        #print("fld_isok:", self.fld_isok)
        if not self.fld_isok or self.status == "failed":
            try:
                stderr_text = str2html(self.stderr_read())
                stdout_text = str2html(tail_file(self.stdout_fname, nlast))
                abiabort_text = "No __ABI_MPIABORTFILE__ found"

                if abort_file:
                    with open(abort_file, "rt") as f:
                        abiabort_text = 12*"=" + os.path.basename(abort_file) + 12*"=" + 2*"\n" + str(f.read())

            except Exception as exc:
                s = "Exception while trying to get info from stderr, stdout and __ABI_MPIABORTFILE\n" + str(exc)
                stderr_text, stdout_text, abiabort_text = 3 * (s,)

            # Look for extra info on the error in selected files produced by the code.
            try:
                errinfo_text = str2html(extract_errinfo_from_files(self.workdir))
            except Exception as exc:
                errinfo_text = "Exception while trying to get error info from extra files\n" + str(exc)

        ##################################################
        # Document Name Space that serves as the substitution
        # namespace for instantiating a doc template.
        try:
            username = os.getlogin()
        except:
            username = "No_username"

        DNS = {
            "self": self,
            "page_title": "page_title",
            "user_name": username,
            "hostname": gethostname(),
            "Headings": ['File_to_test', 'Status', 'fld_output', 'fld_options', 'txt_diff', 'html_diff',] ,
            "nlast": nlast,
            "stderr_text": stderr_text,
            "stdout_text": stdout_text,
            "abiabort_text": abiabort_text,
            "errinfo_text": errinfo_text,
            # Functions and modules available in the template.
            "time": time,
            "pj": os.path.join,
            "basename": os.path.basename,
            "str2html": str2html,
            "sec2str": sec2str,
            "args2htmltr": args2htmltr,
            "html_link"  : html_link,
            "status2html": status2html
            }

        header = """
        <html>
         <head><title>$page_title</title></head>
         <body bgcolor="#FFFFFF" text="#000000">
        """

        if self.status in ["skipped", "disabled"]:
            if self.status == "skipped":
                template = str2html(self.skip_msg)
            else:
                template = "This test has been disabled!"
        else:
            template = """
              <hr>
              <h1>Results of test ${self.full_id}</h1>
                 MPI nprocs =  ${self.nprocs},
                 run_etime = ${sec2str(self.run_etime)} s,
                 tot_etime = ${sec2str(self.tot_etime)} s
               <br>
               ${html_link("stdin",  basename(self.stdin_fname))},
               ${html_link("stdout", basename(self.stdout_fname))},
               ${html_link("stderr", basename(self.stderr_fname))}
              <p>
              <table width="100%" border="0" cellspacing="0" cellpadding="2">
                <tr valign="top" align="left">
                <py-open code = "for h in Headings:"> </py-open>
                  <th>${h}</th>
                <py-close/>
                </tr>
                <py-open>for idx, f in enumerate(self.files_to_test):</py-open>
                 <tr valign="top" align="left">
                  <py-line code = "fld_link = html_link(basename(f.fldiff_fname))"/>
                  <py-line code = "txt_diff_link = html_link(basename(f.diff_fname))"/>
                  <py-line code = "html_diff_link = html_link(basename(f.hdiff_fname))"/>
                  <py-line code = "tab_row = args2htmltr(f.name, status2html(f.fld_status), fld_link, f.fld_options, txt_diff_link, html_diff_link)"/>
                  ${tab_row}
                 </tr>
                <py-close/>
              </table>

              <py-open>for idx, f in enumerate(self.files_to_test):</py-open>
                <py-open code="if f.fld_status != 'succeeded':"/>
                <p> ${f.name} ${f.fld_msg} </p>
              <py-close/>

              <py-open code="if self.status == "failed":"/>
                <py-open code="if self.exceptions:"/>
                  <hr><p>
                  <h1>Exceptions raised at run-time:</h1>
                  <py-open code="for idx, e in enumerate(self.exceptions):"/>
                    <p> $idx) ${str2html(str(e))}</p>
                  <py-close/>
                  <br>
                <py-close/>
                <hr><p>
                <h1>Standard Error of test ${self.id}:</h1>
                  ${stderr_text}
                <hr><p>
                <h1>__MPIABORTFILE__ of test ${self.id}:</h1>
                  ${abiabort_text}
                <hr><p>
                <h1>Info extracted from debug files produced by ${self.id}:</h1>
                  ${errinfo_text}
                <hr><p>
                <h1>Standard output of test ${self.id} (last ${nlast} lines):</h1>
                  ${stdout_text}
                <br>
              <py-close/>
              <p>
              <h3>Extra Information</h3>
              <py-line code = "authors = ', '.join([a for a in self.authors])" />
              <p>Authors = ${authors}</p>
              <py-line code = "keys = ', '.join([k for k in self.keywords])" />
              <p>Keywords = ${keys}</p>
              <p>${self.listoftests(abslink=False)}</p>
            """

        footer = """
          <hr>
          Automatically generated by %s on %s. Logged on as %s@%s
          Python version: %s
          <hr>
          </body>
          </html> """ % (_MY_NAME, time.asctime(), username, gethostname(), platform.python_version())

        if "o" in oc: template = header + template
        if "c" in oc: template += footer

        # Set a file-like object to template
        template_stream = StringIO(template)

        # Initialise an xyaptu xcopier, and call xcopy
        xcp = xcopier(DNS, ouf=fh)
        xcp.xcopy(template_stream)

        if close_fh: fh.close()

    def _get_one_backtrace(self):
        return NagBacktrace(self.stderr_readlines())

    def get_backtraces(self):
        return [self._get_one_backtrace()]

#############################################################################################################
# Subclasses needed to handle the different executables
#############################################################################################################


class AbinitTest(BaseTest):
    """
    Class for Abinit tests. Redefine the make_stdin method of BaseTest
    """
    def make_stdin(self):
        t_stdin = StringIO()

        inp_fname = self.cygwin_path(self.inp_fname)
        # Use the basename instead of the absolute path because the input has been already copied
        # and we might want to change it especially if we are debugging the code
        inp_fname = os.path.basename(inp_fname)
        t_stdin.write(inp_fname + "\n")

        out_fname = self.id + ".out"
        t_stdin.write(out_fname + "\n")

        # Prefix for input-output-temporary files
        if self.input_prefix:
            i_prefix = self.input_prefix
        else:
            i_prefix = self.id + "i"

        # Prefix for input-output-temporary files
        if self.output_prefix:
            o_prefix = self.output_prefix
        else:
            o_prefix = self.id + "o"

        t_prefix = self.id #+ "t"

        t_stdin.writelines([l + "\n" for l in [i_prefix, o_prefix, t_prefix]])

        # Path to the pseudopotential files.
        # 1) pp files are searched in pspd_dir first then in workdir.
        psp_paths = [os.path.join(self.abenv.psps_dir, pname) for pname in self.psp_files]

        for idx, psp in enumerate(psp_paths):
            if not os.path.isfile(psp):
                pname = os.path.join(self.workdir, os.path.basename(psp))
                if os.path.isfile(pname):
                    # Use local pseudo.
                    psp_paths[idx] = pname
                else:
                    err_msg = "Cannot find pp file %s, neither in Psps_for_tests nor in self.workdir" % pname
                    self.exceptions.append(self.Error(err_msg))

        psp_paths = [self.cygwin_path(p) for p in psp_paths] # Cygwin

        t_stdin.writelines([p + "\n" for p in psp_paths])

        return t_stdin.getvalue()


class AnaddbTest(BaseTest):
    """
    Class for Anaddb tests. Redefine the make_stdin method of BaseTest
    """
    def make_stdin(self):
        t_stdin = StringIO()

        inp_fname = self.cygwin_path(self.inp_fname)  # cygwin
        t_stdin.write( inp_fname + "\n")              # 1) formatted input file
        t_stdin.write( self.id + ".out" + "\n")       # 2) formatted output file e.g. t13.out

        iddb_fname = self.id + ".ddb.in"
        if self.input_ddb:
            iddb_fname = os.path.join(self.workdir, self.input_ddb)  # Use output DDB of a previous run.

            if not os.path.isfile(iddb_fname):
                self.exceptions.append(self.Error("%s no such DDB file: " % iddb_fname))

            iddb_fname = self.cygwin_path(iddb_fname)   # cygwin

        t_stdin.write( iddb_fname + "\n")         # 3) input derivative database e.g. t13.ddb.in
        t_stdin.write( self.id + ".md" + "\n")    # 4) output molecular dynamics e.g. t13.md

        input_gkk = self.id + ".gkk"
        if self.input_gkk:
            input_gkk = os.path.join(self.workdir, self.input_gkk) # Use output GKK of a previous run.
            if not os.path.isfile(input_gkk):
                self.exceptions.append(self.Error("%s no such GKK file: " % input_gkk) )

            input_gkk = self.cygwin_path(input_gkk)    # cygwin

        t_stdin.write(input_gkk + "\n")         # 5) input elphon matrix elements  (GKK file) :
        t_stdin.write(self.id + "\n")           # 6) base name for elphon output files e.g. t13

        input_ddk = self.id + ".ddk"
        if not os.path.isfile(input_ddk): # Try in input directory:
            input_ddk = os.path.join(self.inp_dir, input_ddk)
            # FIXME: Someone has to rewrite the treatment of the anaddb files file
            input_ddk = self.cygwin_path(input_ddk)

        t_stdin.write(input_ddk + "\n")   # 7) file containing ddk filenames for elphon/transport :

        return t_stdin.getvalue()

class MultibinitTest(BaseTest):
    """
    Class for Multibinit tests. Redefine the make_stdin method of BaseTest
    """
    def make_stdin(self):
        t_stdin = StringIO()

        inp_fname = self.cygwin_path(self.inp_fname)  # cygwin
        t_stdin.write( inp_fname + "\n")              # 1) formatted input file
        t_stdin.write( self.id + ".out" + "\n")       # 2) formatted output file e.g. t13.out

        if self.input_ddb:
            iddb_fname = os.path.join(self.inp_dir,self.input_ddb)
            if not os.path.isfile(iddb_fname):
                self.exceptions.append(self.Error("%s no such DDB file: " % iddb_fname))
            iddb_fname = self.cygwin_path(iddb_fname)   # cygwin
            t_stdin.write(iddb_fname + "\n")         # 3) input derivative database e.g. ddb.in
        else:
            if self.system_xml:
                sys_xml_fname =  os.path.join(self.inp_dir,self.system_xml)
                if not os.path.isfile(sys_xml_fname):
                    self.exceptions.append(self.Error("%s no such xml file: " % sys_xml_fname))
                sys_xml_fname = self.cygwin_path(sys_xml_fname)
                t_stdin.write(sys_xml_fname + "\n") # 3) input for system.xml XML
            else:
                self.exceptions.append(self.Error("%s no file avail for the system"))

        if self.coeff_xml:
            coeffxml_fname =  os.path.join(self.inp_dir,self.coeff_xml)
            if not os.path.isfile(coeffxml_fname):
                self.exceptions.append(self.Error("%s no such xml file for coeffs: " % coeffxml_fname))

            coeffxml_fname = self.cygwin_path(coeffxml_fname)
            t_stdin.write(coeffxml_fname + "\n") # 4) input for coefficients
        else:
            coeffxml_fname = "no"
            t_stdin.write(coeffxml_fname + "\n")

        if self.md_hist:
            md_hist_fname =  os.path.join(self.inp_dir,self.md_hist)
            if not os.path.isfile(md_hist_fname):
                self.exceptions.append(self.Error("%s no such xml file for coeffs: " % md_hist_fname))

            md_hist_fname = self.cygwin_path(md_hist_fname)
            t_stdin.write(md_hist_fname + "\n") # 5) input for coefficients
        else:
            md_hist_fname = "no"
            t_stdin.write(md_hist_fname + "\n")

        return t_stdin.getvalue()


class AimTest(BaseTest):
    """
    Class for Aim tests. Redefine the make_stdin method of BaseTest
    """
    def make_stdin(self):
        t_stdin = StringIO()

        inp_fname = self.cygwin_path(self.inp_fname)
        t_stdin.write(inp_fname + "\n")         # formatted input file e.g. .../Input/t57.in

        iden_fname = self.id + "i_DEN"
        t_stdin.write(iden_fname + "\n")        # input density  e.g. t57i_DEN
        t_stdin.write(self.id + "\n")           # t57

        # Path to the pseudopotential files.
        psp_paths = [os.path.join(self.abenv.psps_dir, pname) for pname in self.psp_files]
        psp_paths = [self.cygwin_path(p) for p in psp_paths] # Cygwin

        t_stdin.writelines([p + "\n" for p in psp_paths])

        return t_stdin.getvalue()


class ConductiTest(BaseTest):
    """
    Class for Conducti tests. Redefine the make_stdin method of BaseTest
    """
    def make_stdin(self):
        t_stdin = StringIO()

        inp_fname = self.cygwin_path(self.inp_fname)
        t_stdin.write(inp_fname + "\n")  # formatted input file e.g. .../Input/t57.in
        t_stdin.write(self.id + "\n")    # will be used as the prefix of the log file names e.g. t57

        return t_stdin.getvalue()


class OpticTest(BaseTest):
    """
    Class for Optic tests. Redefine the make_stdin method of BaseTest
    """
    def make_stdin(self):
        t_stdin = StringIO()

        inp_fname = self.cygwin_path(self.inp_fname)
        t_stdin.write( inp_fname + "\n")   # optic input file e.g. .../Input/t57.in
        t_stdin.write(self.id + ".out\n")  # Output. e.g t57.out
        t_stdin.write(self.id + "\n")      # Used as suffix to diff and prefix to log file names,
                                           # and also for roots for temporaries

        return t_stdin.getvalue()


class Band2epsTest(BaseTest):
    """How to waste lines of code just to test a F90 code that can be implemented with a few python commands!"""
    def make_stdin(self):
        t_stdin = StringIO()

        inp_fname = self.cygwin_path(self.inp_fname)
        t_stdin.write( inp_fname + "\n")        # input file e.g. .../Input/t51.in
        t_stdin.write( self.id + ".out.eps\n")  # output file e.g. t51.out.eps

        inp_freq = os.path.join(self.inp_dir, self.id + ".in_freq")
        inp_freq = self.cygwin_path(inp_freq)
        t_stdin.write(inp_freq + "\n")          # input freq file e.g Input/t51.in_freq

        inp_displ = os.path.join(self.inp_dir, self.id + ".in_displ")
        if not os.path.isfile(inp_displ):
            inp_displ = "no"
        else:
            inp_displ = self.cygwin_path(inp_displ)
        t_stdin.write(inp_displ + "\n")         # input displ file e.g Input/t51.in_displ

        return t_stdin.getvalue()


class AtompawTest(BaseTest):
    """
    Class for Atompaw tests. Redefine the  methods clean_workdir and bin_path provided by BaseTest
    """
    def clean_workdir(self, other_test_files=None):
        """Keep all atompaw output files."""

    @property
    def bin_path(self):
        """atompaw is not located in src/98_main"""
        return self.build_env.path_of_bin("atompaw")


def exec2class(exec_name):
    """
    Return the test class associated to the executable. Default is BaseTest.
    """
    return {
        "abinit": AbinitTest,
        "anaddb": AnaddbTest,
        "aim": AimTest,
        "conducti": ConductiTest,
        "atompaw": AtompawTest,
        "band2eps": Band2epsTest,
        "optic": OpticTest,
        "multibinit": MultibinitTest,
    }.get(exec_name, BaseTest)


class ChainOfTests(object):
    """
    A list of tests that should be executed together due to inter-dependencies.
    It provides the same interface as the one given by BaseTest
    """
    Error = BaseTestError

    def __init__(self, tests):
        self.tests = tuple([t for t in tests])

        self.inp_dir = tests[0].inp_dir
        self.suite_name = tests[0].suite_name
        #
        # Consistency check.
        for t in tests:
            if self.inp_dir != t.inp_dir or self.suite_name != t.suite_name:
                raise self.Error("All tests should be located in the same directory")

        all_keys = [t.keywords for t in self.tests]
        self.keywords = set()
        for ks in all_keys:
            self.keywords = self.keywords.union(ks)

        all_cpp_vars = [t.need_cpp_vars  for t in self.tests]
        self.need_cpp_vars = set()
        for vs in all_cpp_vars:
            self.need_cpp_vars = self.need_cpp_vars.union(vs)

        self._files_to_keep = []

    def __len__(self):
        return len(self.tests)

    def __str__(self):
        return "\n".join([ str(t) for t in self ])

    def __iter__(self):
        for t in self.tests: yield t

    def info_on_chain(self):
        attr_names = ["extra_inputs", "pre_commands", "post_commands"]
        string = "Info on chain: %s\n" % self.full_id

        nlinks = 0
        for test in self:
            string += test.full_id + "executable " + test.executable + ":\n"
            for (attr, value) in test.__dict__.items():
                if (value and (attr in attr_names or
                    attr.startswith("input_") or attr.startswith("output_"))):
                    string += "  %s = %s\n" % (attr, value)
                    nlinks += 1

        return string, nlinks

    # A lot of boilerplate code!
    # See the doc strings of BaseTest
    @property
    def id(self):
        return "-".join([test.id for test in self])

    @property
    def full_id(self):
        return "["+self.suite_name+"]["+self.id+"]"

    @property
    def max_nprocs(self):
        return max([test.max_nprocs for test in self])

    @property
    def _executed(self):
        return all([test._executed for test in self])

    @property
    def ref_dir(self):
        ref_dirs = [test.ref_dir for test in self]
        assert all([dir == ref_dirs[0] for dir in ref_dirs])
        return ref_dirs[0]

    def listoftests(self, width=100, html=True, abslink=True):
        string = ""
        if not html:
            string += "\n".join( [test.listoftests(width, html, abslink) for test in self] )
            string = self.full_id + ":\n" + string
        else:
            string += "<br>".join( [test.listoftests(width, html, abslink) for test in self] )
            string = "Test Chain " + self.full_id + ":<br>" + string
        return string

    @property
    def files_to_test(self):
        files = []
        for test in self: files.extend(test.files_to_test)
        return files

    @property
    def extra_inputs(self):
        extra_inputs = []
        for test in self: extra_inputs.extend(test.extra_inputs)
        return extra_inputs

    @property
    def inputs_used(self):
        inputs = []
        for test in self: inputs.extend(test.inputs_used)
        return inputs

    @property
    def run_etime(self):
        return sum([test.run_etime for test in self])

    @property
    def tot_etime(self):
        return sum([test.tot_etime for test in self])

    @property
    def isok(self):
        return all([test.isok for test in self])

    @property
    def exceptions(self):
        excs = []
        for test in self:
            excs.extend(test.exceptions)

        return excs

    @property
    def status(self):
        _stats = [test._status for test in self]
        if "disabled" in _stats or "skipped" in _stats:
            if any([s != _stats[0] for s in _stats]):
                #print(self)
                #print("WARNING, expecting all(s == _stats[0] but got\n %s" % str(_stats))
                return "failed"
            return _stats[0]

        all_fldstats = [f.fld_status for f in self.files_to_test]
        if "failed" in all_fldstats: return "failed"
        if "passed" in all_fldstats: return "passed"

        if any([s != "succeeded" for s in all_fldstats]):
            print(self)
            print("WARNING, expecting all(s == 'succeeded' but got\n %s" % str(all_fldstats))

            return "failed"

        return "succeeded"

    def keep_files(self, files):
        if is_string(files):
            self._files_to_keep.append(files)
        else:
            self._files_to_keep.extend(files)

    @property
    def files_to_keep(self):
        # The files produced by the individual tests.
        files_of_tests = []
        for test in self:
            files_of_tests.extend(test.files_to_keep)

        # Add the files produced by self.
        self._files_to_keep += files_of_tests
        return self._files_to_keep

    def cpkl_dump(self, protocol=-1):
        self.cpkl_fname = os.path.join(self.workdir, self.id + ".cpkl")
        with open(self.cpkl_fname, "wb") as fh:
            pickle.dump(self, fh, protocol=protocol)
            self.files_to_keep.append(self.cpkl_fname)

    def has_keywords(self, keywords, mode="any"):
        if mode == "all":
            return set(keywords).issubset(self.keywords)
        elif mode == "any":
            return set(keywords).intersection(self.keywords)
        else:
            raise ValueError("wrong mode %s" % mode)

    def has_variables(self, ivars):
        for test in self:
            matches = test.has_variables(ivars)
            if matches:
                return matches

        return []

    def edit_input(self, editor=None):
        if editor is None: editor = Editor()

        for test in self:
            try:
                test.edit_input(editor=editor)
            except:
                raise

    @property
    def _authors_snames(self):
        snames = set()
        for test in self:
            snames = snames.union(test._authors_snames)
        return snames

    def has_authors(self, authors, mode="any"):
        #return set(authors).issubset(self._authors_snames)
        if mode == "all":
            return set(authors).issubset(self._authors_snames)
        elif mode == "any":
            return set(authors).intersection(self._authors_snames)
        else:
            raise ValueError("wrong mode %s" % mode)

    def write_html_report(self):
        html_report = os.path.join(self.workdir, "test_report.html")
        with open(html_report, "w") as fh:
            for idx, test in enumerate(self):
                oc = ""
                if idx == 0: oc += "o"
                if idx == (len(self)-1): oc += "c"
                test.write_html_report(fh=fh, oc=oc)

    def run(self, build_env, runner, workdir, nprocs=1, **kwargs):

        workdir = os.path.abspath(workdir)
        if not os.path.exists(workdir): os.mkdir(workdir)
        self.workdir = workdir

        for test in self:
            test.run(build_env, runner, workdir=self.workdir, nprocs=nprocs, **kwargs)

    def clean_workdir(self, other_test_files=None):
        for test in self:
            test.clean_workdir(other_test_files=self.files_to_keep)

    def patch(self, patcher=None):
        for test in self:
            test.patch(patcher)

    def get_backtraces(self):
        return [test._get_one_backtrace() for test in self]


class AbinitTestSuite(object):
    """
    List of BaseTest instances. Provide methods to:

    1) select subset of tests according to keywords, authors, numbers
    2) run tests in parallel with python threads
    3) analyze the final results
    """
    def __init__(self, abenv, inp_files=None, test_list=None, keywords=None, need_cpp_vars=None):

        # Check arguments.
        args = [inp_files, test_list]
        no_of_notnone = [arg is not None for arg in args].count(True)
        if no_of_notnone != 1:
            raise ValueError("Wrong args: " + str(args))

        self._executed = False
        self.abenv = abenv
        self.exceptions = []

        if inp_files is not None:
            self.tests = make_abitests_from_inputs(
                inp_files, abenv,
                keywords=keywords, need_cpp_vars=need_cpp_vars)

        elif test_list is not None:
            assert keywords is None
            assert need_cpp_vars is None
            self.tests = tuple(test_list)

        else:
            raise ValueError("Either inp_files or test_list must be specified!")

    def __str__(self):
        return "\n".join([str(t) for t in self.tests])

    def __add__(self, other):
        test_list = [t for t in self] + [t for t in other]
        return self.__class__(self.abenv, test_list=test_list)

    def __len__(self):
        return len(self.tests)

    def __iter__(self):
        for t in self.tests:
            yield t

    def __getitem__(self, key):   # FIXME: this won't work for tutorial, paral and other test suites.
        """Called by self[key]."""
        if isinstance(key, slice):
            return self.__getslice(key)
        else:
            raise NotImplementedError("__getitem__ expects a slice instance")

    def __getslice(self, slice):
        start = slice.start
        if start is None: start = 0
        stop = slice.stop
        if stop is None: stop = 10000 # Not very elegant, but cannot use len(self) since indices are not contiguous
        assert slice.step is None     # Slices with steps (e.g. [1:4:2]) are not supported.

        # Rules for the test id:
        # Simple case: t01, tgw1_1
        # test chain (no MPI): t81-t82-t83-t84, tudet_1-tudet_2-tudet_3
        # multi-parallel tests:  t74_MPI2, t51_MPI1-t52_MPI1-t53_MPI1, tdfpt_01_MPI2 ...

        test_list = []
        for test in self:
            #print("ID",test.id)
            # extract the ID of the first test (if test_chain)
            tokens = test.id.split("-")
            assert tokens[0][0] == "t"  # Assume first character is "t"
            num = tokens[0][1:]

            if "_MPI" in test.id:
                # Handle multi-parallel tests.
                #print(test.id)
                idx = test.id.find("_MPI")
                tok = test.id[1:idx]
                #print(tok)
                idx = tok.rfind("_")
                if idx != -1:
                    # Handle tdfpt_01_MPI2 ...
                    # FIXME: this will fail if _OMP2_MPI2
                    tok = tok[idx+1:]
                try:
                    num = int(tok)
                except ValueError:
                    raise ValueError("Cannot convert %s to integer" % tok)

            else:
                # Simple case or test_chain
                idx = num.rfind("_")
                if idx != -1:
                    num = int(num[idx+1:])

            num = int(num)
            if num in range(start, stop):
                #print "got", test.id
                test_list.append(test)

        return self.__class__(self.abenv, test_list=test_list)

    @property
    def full_length(self):
        one = lambda : 1
        return sum([getattr(test, "__len__", one)() for test in self])

    @property
    def run_etime(self):
        """Total elapsed time i.e. the wall-time spent in the sub-processes e.g. abinit)"""
        assert self._executed
        return sum([test.run_etime for test in self])

    @property
    def keywords(self):
        keys = []
        for test in self: keys.extend(test.keywords)
        return set(keys)

    def has_keywords(self, keywords):
        return set(keywords).issubset(self.keywords)

    @property
    def need_cpp_vars(self):
        vars = []
        for test in self: vars.extend(test.need_cpp_vars)
        return set(vars)

    def on_refslave(self):
        """True if we are running on a reference slave e.g. testf."""
        try:
            return self._on_ref_slave
        except AttributeError:
            return False

    def set_on_refslave(self, value=True):
        """Attribute setter"""
        self._on_ref_slave = bool(value)

    def all_exceptions(self):
        """Return my exceptions + test exceptions."""
        all_excs = self.exceptions
        for test in self:
            all_excs.extend(test.exceptions)
        return all_excs

    def cpkl_dump(self, protocol=-1):
        self.cpkl_fname = os.path.join(self.workdir, "test_suite.cpkl")
        with open(self.cpkl_fname, "wb") as fh:
            pickle.dump(self, fh, protocol=protocol)

    def _tests_with_status(self, status):
        assert status in BaseTest._possible_status
        #assert self._executed
        return [test for test in self if test.status == status]

    def succeeded_tests(self): return self._tests_with_status("succeeded")
    def passed_tests(self):    return self._tests_with_status("passed")
    def failed_tests(self):    return self._tests_with_status("failed")
    def skipped_tests(self):   return self._tests_with_status("skipped")
    def disabled_tests(self):  return self._tests_with_status("disabled")

    @property
    def targz_fname(self):
        """
        Location of the tarball file with the results in HTML format
        Returns None if the tarball has not been created.
        """
        try:
            return self._targz_fname
        except:
            return None

    def create_targz_results(self):
        """Create the tarball file results.tar.gz in the working directory."""
        assert self._executed
        exclude_exts = [".cpkl", ".py", "pyc", ]

        self._targz_fname = None
        ofname = os.path.join(self.workdir,"results.tar.gz")

        # The most delicate part here is the treatment of the exceptions
        # since the test might not have produced the reference files
        # we want to copy to the server. If something goes wrong, we simply
        # register the exception and we continue the execution.
        try:
            targz = tarfile.open(ofname, "w:gz")

            for test in self:

                # Don't try to collect files for tests that are disabled or skipped.
                if test.status in ["disabled", "skipped"]:
                    continue

                files = set(test.files_to_keep)
                save_files = [f for f in files if not has_exts(f, exclude_exts)]
                #print(save_files)

                # Store stdout files only if the test failed.
                important_status = ["failed",]

                # Special treatment for reference machines
                if self.on_refslave:
                    important_status = ["passed", "failed",]

                if test.status not in important_status:
                    if isinstance(test, ChainOfTests):
                        for t in test:
                            #print "Removing Test Chain", t.stdout_fname
                            save_files.remove(t.stdout_fname)
                    else:
                        #print "Removing", test.stdout_fname
                        save_files.remove(test.stdout_fname)

                for p in save_files:
                    #if not os.path.exists(os.path.join(self.workdir, p)): continue
                    # /foo/bar/suite_workdir/test_workdir/file --> test_workdir/t01/file
                    rpath = os.path.relpath(p, start=self.workdir)
                    #arcname = str(rpath.encode("ascii", "ignore"))
                    arcname = str(rpath)
                    try:
                        #print("targz.add: adding:", p," with arcname ", arcname)
                        #print(type(p), type(arcname))
                        targz.add(p, arcname=arcname)
                    except:
                        # Handle the case in which the output file has not been produced.
                        exc = sys.exc_info()[1]
                        warnings.warn("exception while adding %s to tarball:\n%s" % (p, exc))
                        self.exceptions.append(exc)

            targz.close()

            # Save the name of the tarball file.
            self._targz_fname = ofname

        except:
            exc = sys.exc_info()[1]
            warnings.warn("exception while creating tarball file: %s" % str(exc))
            self.exceptions.append(exc)

    def sanity_check(self):
        all_full_ids = [test.full_id for test in self]
        if len(all_full_ids) != len(set(all_full_ids)):
            raise ValueError("Cannot have more than two tests with the same full_id")

    def run_tests(self, build_env, workdir, runner, nprocs=1, nthreads=1, runmode="static", **kwargs):
        """
        Execute the list of tests (main entry point for client code)

        Args:
            build_env:
                `BuildEnv` instance with info on the build environment.
            workdir:
                Working directory (string)
            runner:
                `JobRunner` instance
            nprocs:
                number of MPI processes to use for a single test.
            nthreads:
                number of OpenMP threads for tests
        """
        self.sanity_check()

        if len(self) == 0:
            warnings.warn("len(self) == 0")
            return

        workdir = os.path.abspath(workdir)
        if not os.path.exists(workdir): os.mkdir(workdir)

        self.workdir = workdir

        # Acquire the lock file.
        self.lock = FileLock(os.path.join(self.workdir,"__run_tests_lock__"), timeout=3)

        try:
            self.lock.acquire()
        except self.lock.Error:
            msg = ("Timeout occured while trying to acquire lock in:\n\t%s\n"
                   "Perhaps a previous run did not exit cleanly or another process is running in the same directory.\n"
                   "If you are sure no other process is in execution, remove the directory with `rm -rf` and rerun.\n" % self.workdir)
            cprint(msg, "red")
            return

        # Remove all stale files present in workdir (except the lock!)
        rmrf(self.workdir, exclude_paths=self.lock.lockfile)

        self.nprocs = nprocs
        self.nthreads = nthreads

        def run_and_check_test(test):
            """Helper function to execute the test. Must be thread-safe."""

            testdir = os.path.abspath(os.path.join(self.workdir, test.suite_name + "_" + test.id))

            # Run the test
            test.run(build_env, runner, testdir, nprocs=nprocs, runmode=runmode, **kwargs)

            # Write HTML summary
            test.write_html_report()

            # Dump the object with pickle.
            #test.cpkl_dump()

            # Remove useless files in workdir.
            test.clean_workdir()

        ##############################
        # And now let's run the tests
        ##############################
        start_time = time.time()

        if nthreads == 1:
            logger.info("Sequential version")
            for test in self:
                run_and_check_test(test)

        elif nthreads > 1:
            logger.info("Threaded version with nthreads = %s" % nthreads)
            from threading import Thread
            # Internal replacement that provides task_done, join_with_timeout (py2.4 compliant)
            from pymods.myqueue import QueueWithTimeout

            def worker():
                while True:
                    test = q.get()
                    run_and_check_test(test)
                    q.task_done()

            q = QueueWithTimeout()
            for i in range(nthreads):
                t = Thread(target=worker)
                t.setDaemon(True)
                t.start()

            for test in self: q.put(test)

            # Block until all tasks are done. Raise QueueTimeoutError after timeout seconds.
            timeout_1test = float(runner.timebomb.timeout)
            if timeout_1test <= 0.1: timeout_1test = 240.

            queue_timeout = 1.3 * timeout_1test * self.full_length / float(nthreads)
            q.join_with_timeout(queue_timeout)

        # Run completed.
        self._executed = True

        # Collect HTML files in a tarball
        self.create_targz_results()

        nsucc = len(self.succeeded_tests())
        npass = len(self.passed_tests())
        nfail = len(self.failed_tests())
        nskip = len(self.skipped_tests())
        ndisa = len(self.disabled_tests())

        self.tot_etime = time.time() - start_time

        # Print summary table.
        stats_suite = {}
        for test in self:
            if test.suite_name not in stats_suite:
                d = dict.fromkeys(BaseTest._possible_status, 0)
                d["run_etime"] = 0.0
                d["tot_etime"] = 0.0
                stats_suite[test.suite_name] = d

        for test in self:
            stats_suite[test.suite_name][test.status] += 1
            stats_suite[test.suite_name]["run_etime"] += test.run_etime
            stats_suite[test.suite_name]["tot_etime"] += test.tot_etime

        suite_names = sorted(stats_suite.keys())

        times = ["run_etime", "tot_etime"]

        table = [["Suite"] + BaseTest._possible_status + times]
        for suite_name in suite_names:
            stats = stats_suite[suite_name]
            row = [suite_name] + [str(stats[s]) for s in BaseTest._possible_status] + ["%.2f" % stats[s] for s in times]
            table.append(row)

        print("")
        pprint_table(table)
        print("")

        executed = [t for t in self if t.status != "skipped"]
        if executed:
            mean_etime = sum([test.run_etime for test in executed]) / len(executed)
            dev_etime = (sum([(test.run_etime - mean_etime)**2 for test in executed]) / len(executed))**0.5

            cprint("Completed in %.2f [s]. Average time for test=%.2f [s], stdev=%.2f [s]" % (
                  self.tot_etime, mean_etime, dev_etime), "yellow")

            msg = "Summary: failed=%s, succeeded=%s, passed=%s, skipped=%s, disabled=%s" % (
                  nfail, nsucc, npass, nskip, ndisa)

            if nfail:
              cprint(msg, "red", attrs=['underline'])
            else:
              cprint(msg, "green")

            # Print outliers
            if False and dev_etime > 0.0:
                for test in self:
                    if abs(test.run_etime) > 0.0 and abs(test.run_etime - mean_etime) > 2 * dev_etime:
                        print("%s has run_etime %.2f s" % (test.full_id, test.run_etime))

        with open(os.path.join(self.workdir, "results.txt"), "w") as fh:
            pprint_table(table, out=fh)

        try:
            username = os.getlogin()
        except:
            username = "No_username"

        # Create the HTML index.
        DNS = {
            "self": self,
            "runner": runner,
            "user_name": username,
            "hostname": gethostname(),
            "test_headings": ['ID', 'Status', 'run_etime (s)', 'tot_etime (s)'],
            "suite_headings": ['failed', 'passed', 'succeeded', 'skipped', 'disabled'],
            # Functions and modules available in the template.
            "time": time,
            "pj": os.path.join,
            "basename": os.path.basename,
            "str2html": str2html,
            "sec2str": sec2str,
            "args2htmltr": args2htmltr,
            "html_link": html_link,
            "status2html": status2html,
             }

        fname = os.path.join(self.workdir, "suite_report.html")
        fh = open(fname, "w")

        header = """
          <html>
          <head><title>Suite Summary</title></head>
          <body bgcolor="#FFFFFF" text="#000000">
           <hr>
           <h1>Suite Summary</h1>
            <table width="100%" border="0" cellspacing="0" cellpadding="2">
            <tr valign="top" align="left">
             <py-open code = "for h in suite_headings:"> </py-open>
             <th>${status2html(h)}</th>
            <py-close/>
            </tr>
            <tr valign="top" align="left">
            <py-open code = "for h in suite_headings:"> </py-open>
             <td> ${len(self._tests_with_status(h))} </td>
            <py-close/>
            </tr>
            </table>
            <p>
            tot_etime = ${sec2str(self.tot_etime)} <br>
            run_etime = ${sec2str(self.run_etime)} <br>
            no_pythreads = ${self.nthreads} <br>
            no_MPI = ${self.nprocs} <br>
            ${str2html(str(runner))}
           <hr>
        """

        table = """
           <p>
           <h1>Test Results</h1>
           <table width="100%" border="0" cellspacing="0" cellpadding="2">
            <tr valign="top" align="left">
             <py-open code = "for h in test_headings:"> </py-open>
              <th>$h</th>
             <py-close/>
            </tr>
        """

        for status in BaseTest._possible_status:
            table += self._pyhtml_table_section(status)

        table += "</table>"

        footer = """
          <hr>
          <h1>Suite Info</h1>
            <py-line code = "keys = ', '.join(self.keywords)" />
            <p>Keywords = ${keys}</p>
            <py-line code = "cpp_vars = ', '.join(self.need_cpp_vars)"/>
            <p>Required CPP variables = ${cpp_vars}</p>
          <hr>
            Automatically generated by %s on %s. Logged on as %s@%s
          <hr>
          </body>
          </html> """ % (_MY_NAME, time.asctime(), username, gethostname() )

        template = header + table + footer

        template_stream = StringIO(template)

        # Initialise an xyaptu xcopier, and call xcopy
        xcp = xcopier(DNS, ouf=fh)
        xcp.xcopy(template_stream)
        fh.close()

        # Release the lock.
        self.lock.release()

        return Results(self)

    @staticmethod
    def _pyhtml_table_section(status):
        # ['ID', 'Status', 'run_etime', 'tot_etime'],
        string = """
           <py-open code="for test in self.%s_tests():"/>
            <py-line code = "report_link = pj(basename(test.workdir),'test_report.html') " />
            <tr valign="top" align="left">
             <td> ${html_link(test.full_id, report_link)}</td>
             <td> ${status2html(test.status)} </td>
             <td> ${sec2str(test.run_etime)} </td>
             <td> ${sec2str(test.tot_etime)} </td>
            </tr>
           <py-close/>
           """ % status
        return string

    def patch(self, patcher=None):
        """
        Patch the output files of the test with the specified patcher.
        A default patcher is provided if patcher is None (use $PATCHER shell variable)
        """
        for test in self:
            test.patch(patcher)

    def select_tests(self, with_keys=None, exclude_keys=None, with_authors=None, exclude_authors=None,
                     ivars=None, mode="any"):
        """
        Extract the subset of tests matching the given conditions.

        Returns:
            `AbinitTestSuite` instance
        """
        test_list = [test for test in self]

        if with_keys:
            test_list = [test for test in test_list if test.has_keywords(with_keys, mode=mode)]

        if exclude_keys:
            test_list = [test for test in test_list if not test.has_keywords(exclude_keys, mode=mode)]

        if with_authors:
            test_list = [test for test in test_list if test.has_authors(with_authors, mode=mode)]

        if exclude_authors:
            test_list = [test for test in test_list if not test.has_authors(exclude_authors, mode=mode)]

        if ivars:
            test_list = [test for test in test_list if test.has_variables(ivars)]

        return AbinitTestSuite(self.abenv, test_list=test_list)

    def make_listoftests(self, width=100, html=True):
        """Create the ListOfTests files."""
        if not html:
            return "\n\n".join([test.listoftests(width, html) for test in self])
        else:
            header = """
             <html>
             <head><title>"LIST OF TESTS" FILE</title></head>
             <body bgcolor="#FFFFFF" text="#000000">
             <!-- Automatically generated by %s on %s. ****DO NOT EDIT**** -->""" % (_MY_NAME, time.asctime())

            body = "<hr>".join([test.listoftests(width, html) for test in self])

            footer = """
              <hr>
               Automatically generated by %s on %s.
              <hr>
              </body>
              </html>""" % (_MY_NAME, time.asctime())

            return header + body + footer


class Results(object):
    """Stores the final results."""
    def __init__(self, test_suite):
        #assert test_suite._executed
        self.test_suite = test_suite
        self.failed_tests = test_suite.failed_tests()
        self.passed_tests = test_suite.passed_tests()
        self.succeeded_tests = test_suite.succeeded_tests()
        self.skipped_tests = test_suite.skipped_tests()
        self.disabled_tests = test_suite.disabled_tests()
        self.targz_fname = test_suite.targz_fname

    @lazy__str__
    def __str__(self): pass

    def tests_with_status(self, status):
        return {
            "succeeded": self.succeeded_tests,
            "passed": self.passed_tests,
            "failed": self.failed_tests,
            "disabled": self.disabled_tests,
            "skipped": self.skipped_tests,
            "all": [test for test in self.test_suite]
            }[status]

    @property
    def nfailed(self):
        """Number of failures"""
        return len(self.failed_tests)

    @property
    def npassed(self):
        """Number of tests marked as passed."""
        return len(self.passed_tests)

    @property
    def nexecuted(self):
        """Number of tests executed."""
        n = 0
        for test in self.test_suite:
            if isinstance(test, ChainOfTests):
                n += len([t for t in test if t.executed])
            else:
                if test.executed:
                    n += 1
        return n

    def outref_files(self, status):
        """
        Return (out_files, ref_files)
        where out and ref are list with the output files and the reference
        files of the tests with the given status.
        """
        out_files, ref_files = [], []
        for test in (self.tests_with_status(status)):
            for f in test.files_to_test:
                #if status != "all" and f.status != status: continue

                out_files.append(os.path.join(test.workdir, f.name))
                ref_fname = os.path.join(test.ref_dir, f.name)
                # FIXME Hack due to the ambiguity stdout, out!
                if not os.path.exists(ref_fname) and ref_fname.endswith(".stdout"):
                    ref_fname = ref_fname[:-7] + ".out"
                ref_files.append(ref_fname)

        return out_files, ref_files

    def in_files(self, status):
        """List with the input files of the tests with the given status."""
        in_files = []
        for test in (self.tests_with_status(status)):
            if isinstance(test, ChainOfTests):
                in_files.extend([t.inp_fname for t in test])
            else:
                in_files.append(test.inp_fname)

        return in_files

    def patch_refs(self, status="failed"):
        """Patch the reference files of the tests with the specified status."""
        out_files, ref_files = self.outref_files(status=status)
        #for r, o in zip(out_files, ref_files): print("reference: %s, output %s" % (r, o))

        return Patcher().patch_files(out_files, ref_files)

    def edit_inputs(self, status="failed"):
        """Edit the input files of the tests with the specified status."""
        in_files = self.in_files(status=status)
        #for r, o in zip(out_files, ref_files):
        #    print("reference: %s, output %s" % (r, o))

        return Editor().edit_files(in_files)

    #def inspect_stdouts(self):
    #  out_files, ref_files = self.outref_files()
    #  return Editor().edit_files(in_files)
    #def inspect_diffs(self):

    def inspect_stderrs(self, status="failed"):
        """Open the stderr of the tests with the give status in `Editor`."""
        return Editor().edit_files(self.stderr_files(status))

    def stderr_files(self, status="failed"):
        """List of non-empty error files of the tests with the specified status."""
        # Loop over the tests, open the stderr to see if it's empty ot not
        # and add it to the list.
        err_files = []
        for test in self.tests_with_status(status):
            if isinstance(test, ChainOfTests):
                es = [t.stderr_fname for t in test if not t.has_empty_stderr]
                if es:
                    err_files.extend(es)
            else:
                if not test.has_empty_stderr: err_files.append(test.stderr_fname)

        return err_files

    def cpkl_dump(self, cpkl_fname, protocol=-1):
        """Save the object in pickle format."""
        with open(cpkl_fname, "wb") as fh:
            pickle.dump(self, fh, protocol=protocol)


if __name__ == "__main__":
    # Automatic documentation of the TEST_INFO options.
    doc_testcnf_format()
