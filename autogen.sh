aclocal -I m4
autoheader
libtoolize -f
autoconf
automake -a -c
./configure LOG4CXX_CFLAGS="`pkg-config --cflags liblog4cxx`" LOG4CXX_LIBS="`pkg-config --libs liblog4cxx`"
