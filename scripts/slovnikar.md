# Gama-local translations

The program ```gama-local``` supports textual output in various
languages. The language is specified using the ```--language``` option
(the default is English). Language codes follow the two-letter ISO 639
standard. Translations are defined in the following files:

    lib/gnu_gama/local/language.[h|cpp]

These files are generated by the program ```scripts/slovnikar.cpp```,
which contains all the essential information in its comments. The best
way to add a new language is by running the script:

    scripts/build-dictionaries.sh

Translations are defined in several XML files from which files
```language.[h|cpp]``` are compiled. Here is a trivial example (file
```xml/lang/en/gamalib.lang```)

	<?xml version="1.0" encoding="utf-8" standalone="yes" ?>
	<entries>

    <e id="T_IE_internal_error"
	   en="internal error" />

	<e id="T_GaMa_from_equals_to"
	   en="station and target are identical" />

	</entries>

 The meaning is (hopefully) self evident. Entries are defined inside
 tags ```<e />``` with compulsory attribute ```id```, which identifies
 the translated text token whose translations/definition is given in
 attributes corresponding to ISO 639 language names, here ```en```. It
 would be possible to use two or more translations within a single
 ```<e />``` but in our project we keep translations separated (they
 are compiled together by an auxiliary program
 ```scripts/slovnikar.cpp```).

To ease the process of language support build we suggest to use a
shell script ```scripts/build-dictionaries.sh``` run from the top
project directory.

	#!/bin/sh
	#
	# builds dictionaris for GNU Gama
	#

	echo
	echo Building dictionaries
	echo

	if pwd | grep /scripts$; then cd ..; fi

	SLOVNIKAR=$(pwd)/scripts/slovnikar

	# source directory for language.{h|cpp}
	cd lib/gnu_gama/local
	FILES=$(find ../../../xml/lang -name *.lang | grep -v 00.lang)

	echo input files: $FILES
	echo
	echo $FILES | $SLOVNIKAR

**NOTE**: We expect the project to be build as ```./autogen.sh & configure &
make``` from the top project direcory.

**NOTE**: Together with files ```language.[h|cpp]``` in
```lib/gnu_gama_local``` also a file ```language.html``` is generated
which gives a review of all translated texts available. If some
translation of a message is missing, English equivalent is used,



## Template translation files


Translation files, used for generating source files
```language.[h|cpp]```, are stored In the subdirectory ```xml/lang```.
Each language has its own subdirectory, named by the ISO 639 language
name, with translation files ```*.lang``` in it.  Names of the file
are composed of major name, dot, language name and ```.lang```
suffix. For example ```gamalib.cz.lang```


	<?xml version="1.0" encoding="utf-8" standalone="yes" ?>
	<entries>

	<e id="T_IE_internal_error"
	   EN="internal error"
	   cz="interní chyba" />

	<e id="T_GaMa_from_equals_to"
	   EN="station and target are identical"
	   cz="stanovisko je totožné s cílem" />

	</entries>

In this example, the attribute ```id``` serves as a unique identifier for a
text unit, while ```cz``` represents the label for its Czech
translation. The XML attribute ```EN``` (written in uppercase) functions as
a comment, containing the English text associated with the given
```id```. This text typically matches the English prototype defined
elsewhere, which, in this case, is located in the file
```xml/lang/en/gamalib.lang```.

	<?xml version="1.0" encoding="utf-8" standalone="yes" ?>
	<entries>

	<e id="T_IE_internal_error"
	   en="internal error" />

	<e id="T_GaMa_from_equals_to"
	   en="station and target are identical" />

	</entries>

It would be tempting to use the prototype English texts, without
introducing somehow clumsy ```id``` attributes, but this would mean to
always uppdate all translation definition when the oritignal English text
needed to be changed. Thus our project relies on a bit clumsy ```id```
attributes with ```EN``` informative comments to help the translators
to prepare their translations.

### Discarding selected id

There is a special attribute ```XX``` for discarding selected
translations of a given ```id```.


	<?xml version="1.0" encoding="utf-8" standalone="yes" ?>
	<entries>

	<!-- discarded id for all languages (attribute XX) -->
	<e id="T_POBS_computation_of_bearing_for_identical_points"
	   XX="Computation of bearing for identical points" />

	<e id="T_POBS_bad_data"
	   en="bad data" />

	<e id="T_POBS_zero_or_negative_distance"
	   en="zero or negative distance" />

	<e id="T_POBS_zero_or_negative_slope_distance"
	   en="zero or negative slope distance" />

	<e id="T_POBS_zero_or_negative_zenith_angle"
	   en="zero or negative zenith angle" />

	</entries>

Attribute ```XX``` with the given ```id``` can be used in any language
file used.