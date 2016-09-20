#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
exe_name=$0
exe_dir=`dirname "$0"`
echo "------------------------------------------"
if [ "x$1" = "x" ]; then
  echo Usage:
  echo    $0 \<deployedMCRroot\> args
else
  echo Setting up environment variables
  MCRROOT="$1"
  echo ---
  LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
	MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64 ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads ; 
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE} ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b/runtime/glnxa64;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b/sys/os/glnxa64;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b/sys/java/jre/glnxa64/jre/lib/amd64/server;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b/lib:/broad/lsf/7.0/linux2.6-glibc2.3-x86_64/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/lsf_drmaa-1.0.3/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.6.5/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/graphviz_2.26.3/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/db_4.7.25/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/scons_1.2.0-python-2.6.5/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/gtk+_2.18.6/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/fontconfig_2.8.0/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/freetype_2.3.11/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/atk_1.28.0/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/pango_1.26.2/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/cairo_1.8.8/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/pixman_0.16.4/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/glib_2.22.4/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/gettext_0.17/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/libiconv_1.13.1/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/biopython_1.53-python-2.6.5/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/django_1.1.1-python-2.6.5/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/matplotlib_0.99.1.1-python-2.6.5/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/scipy_0.7.1-python-2.6.5/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/numpy_1.5.1-python-2.6.5/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/subversion_1.6.5/lib:/broad/software/free/Linux/redhat_5_x86_64/pkgs/sqlite_3.6.19/lib:/broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b/runtime/glnxa64:/broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/matlab_2010b/bin/glnxa64
  XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
  export LD_LIBRARY_PATH;
  export XAPPLRESDIR;
  echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
  shift 1
  args=
  while [ $# -gt 0 ]; do
      token=`echo "$1" | sed 's/ /\\\\ /g'`   # Add blackslash before each blank
      args="${args} ${token}" 
      shift
  done
  eval "${exe_dir}"/annotate_with_conservartion $args
fi
exit

