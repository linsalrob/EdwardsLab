LIBS_PATH=./libs

ifdef DEBUG
	FLAGS=-DDEBUG=1 --debug
else
	FLAGS=-O3
endif

UTILITY_LIB_PATH := $(LIBS_PATH)/utility_lib
STRING_BUF_PATH := $(LIBS_PATH)/string_buffer
BIOINF_LIB_PATH := $(LIBS_PATH)/bioinf
SCORING_PATH := $(LIBS_PATH)/alignment_scoring

# Check mac/linux
UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
	FLAGS := $(FLAGS) -fnested-functions
endif

# Add compile time
FLAGS := $(FLAGS) -DCOMPILE_TIME='"$(shell date)"' -DSCORE_TYPE='int'

all:
	gcc -o needleman_wunsch $(FLAGS) -Wall \
	-I . -I $(UTILITY_LIB_PATH) -I $(STRING_BUF_PATH) \
	-I $(BIOINF_LIB_PATH) -I $(SCORING_PATH) \
	nw_cmdline.c needleman_wunsch.c \
	$(SCORING_PATH)/*.c $(UTILITY_LIB_PATH)/utility_lib.c \
	$(BIOINF_LIB_PATH)/bioinf.c $(STRING_BUF_PATH)/string_buffer.c -lz

clean:
	if test -e needleman_wunsch; then rm needleman_wunsch; fi
	for file in $(wildcard *.dSYM); do rm -r $$file; done
	for file in $(wildcard *.greg); do rm $$file; done
