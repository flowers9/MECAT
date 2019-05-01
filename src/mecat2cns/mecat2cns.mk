ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := mecat2cns
SOURCES  := main.cpp \
	argument.cpp \
	dw.cpp \
	MECAT_AlnGraphBoost.C \
	mecat_correction.cpp \
	options.cpp \
	overlaps_partition.cpp \
	reads_correction_aux.cpp \
	reads_correction_can.cpp \
	reads_correction_m4.cpp \
	packed_db.cpp \

SRC_INCDIRS  := . libboost

TGT_CXXFLAGS := -D _FILE_OFFSET_BITS=64 -pg -g
TGT_LDFLAGS  := -L${TARGET_DIR} -pg -g
TGT_LDLIBS   := -lmecat
TGT_PREREQS  := libmecat.a

SUBMAKEFILES :=
