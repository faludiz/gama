# https://www.gnu.org/software/automake/manual/html_node/Tests.html
#
# Automake can generate code to handle two kinds of test suites. One
# is based on integration with the dejagnu framework. The other (and
# most used) form is based on the use of generic test scripts, and its
# activation is triggered by the definition of the special TESTS
# variable. This second form allows for various degrees of
# sophistication and customization; in particular, it allows for
# concurrent execution of test scripts, use of established test
# protocols such as TAP, and definition of custom test drivers and
# test runners.
#
# In either case, the testsuite is invoked via ‘make check’.

TESTSA = gama-g3-adjustment.sh \
         gama-g3-ellipsoid-xyz2blh.sh \
         gama-g3-ellipsoid-xyz2blh-list.sh

if GNU_GAMA_LOCAL_TEST_XMLLINT
TESTXML = gama-g3-xmllint-xsd.sh
endif

TESTS = $(TESTSA) $(TESTXML)

SUBDIRS = input src

INPUT_FILES = demo-g3-01 demo-g3-02 demo-g3-03 \
              ghilani-gnss-v1 ghilani-gnss-v2 ghilani-gnss-v3

G3_SRC     = $(srcdir)/src
G3_INPUT   = $(srcdir)/input
G3_RESULTS = results/$(PACKAGE_VERSION)
G3_OTHERS  = Makefile $(G3_RESULTS)

do_subst = sed -e 's,[@]top_builddir[@],$(top_builddir),g' \
               -e 's,[@]PACKAGE_VERSION[@],$(PACKAGE_VERSION),g' \
               -e 's,[@]INPUT_FILES[@],$(INPUT_FILES),g' \
               -e 's,[@]G3_INPUT[@],$(G3_INPUT),g' \
               -e 's,[@]G3_RESULTS[@],$(G3_RESULTS),g' \
               -e 's,[@]GAMA_XML[@],$(top_srcdir)/xml,g'

clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -rf $(G3_RESULTS) $(TESTS)

$(G3_RESULTS):
	mkdir -p $(G3_RESULTS)

gama-g3-adjustment.sh: $(G3_SRC)/gama-g3-adjustment.in $(G3_OTHERS)
	@$(do_subst) < $(G3_SRC)/gama-g3-adjustment.in \
	        > gama-g3-adjustment.sh
	@chmod +x gama-g3-adjustment.sh

gama-g3-xmllint-xsd.sh: $(G3_SRC)/gama-g3-xmllint-xsd.in $(G3_OTHERS)
	@$(do_subst)  < $(G3_SRC)/gama-g3-xmllint-xsd.in \
	        > gama-g3-xmllint-xsd.sh
	@chmod +x gama-g3-xmllint-xsd.sh

gama-g3-ellipsoid-xyz2blh.sh: $(G3_SRC)/gama-g3-ellipsoid-xyz2blh.in \
			      $(G3_OTHERS)
	@$(do_subst)  < $(G3_SRC)/gama-g3-ellipsoid-xyz2blh.in \
	        > gama-g3-ellipsoid-xyz2blh.sh
	@chmod +x gama-g3-ellipsoid-xyz2blh.sh

gama-g3-ellipsoid-xyz2blh-list.sh: $(G3_SRC)/gama-g3-ellipsoid-xyz2blh-list.in \
			      $(G3_OTHERS)
	@$(do_subst)  < $(G3_SRC)/gama-g3-ellipsoid-xyz2blh-list.in \
	        > gama-g3-ellipsoid-xyz2blh-list.sh
	@chmod +x gama-g3-ellipsoid-xyz2blh-list.sh
