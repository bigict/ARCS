aclocal
autoheader
case `uname` in
    Darwin*) glibtoolize -f
        ;;
    *) libtoolize -f
esac
autoconf
automake -a -c
