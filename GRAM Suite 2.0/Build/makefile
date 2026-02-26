.PHONY: all clean shared Venus Earth Mars Jupiter Neptune Saturn Titan Uranus GRAM

all:
	$(MAKE) -C googletest all
	$(MAKE) -C GRAMCommon all
ifneq ($(wildcard Venus),)
	$(MAKE) -C Venus all
endif
ifneq ($(wildcard Earth),)
	$(MAKE) -C Earth all
endif
ifneq ($(wildcard Mars),)
	$(MAKE) -C Mars all
endif
ifneq ($(wildcard Jupiter),)
	$(MAKE) -C Jupiter all
endif
ifneq ($(wildcard Neptune),)
	$(MAKE) -C Neptune all
endif
ifneq ($(wildcard Saturn),)
	$(MAKE) -C Saturn all
endif
ifneq ($(wildcard Titan),)
	$(MAKE) -C Titan all
endif
ifneq ($(wildcard Uranus),)
	$(MAKE) -C Uranus all
endif
	$(MAKE) -C GRAM all

GRAM: all

shared:
	$(MAKE) -C googletest all
	$(MAKE) -C GRAMCommon all
ifneq ($(wildcard Venus),)
	$(MAKE) -C Venus all
endif
ifneq ($(wildcard Earth),)
	$(MAKE) -C Earth all
endif
ifneq ($(wildcard Mars),)
	$(MAKE) -C Mars all
endif
ifneq ($(wildcard Jupiter),)
	$(MAKE) -C Jupiter all
endif
ifneq ($(wildcard Neptune),)
	$(MAKE) -C Neptune all
endif
ifneq ($(wildcard Saturn),)
	$(MAKE) -C Saturn all
endif
ifneq ($(wildcard Titan),)
	$(MAKE) -C Titan all
endif
ifneq ($(wildcard Uranus),)
	$(MAKE) -C Uranus all
endif
	$(MAKE) -C GRAM shared
	
Venus:
ifneq ($(wildcard Venus),)
	$(MAKE) -C googletest all
	$(MAKE) -C GRAMCommon all
	$(MAKE) -C Venus all
endif

Earth:
ifneq ($(wildcard Earth),)
	$(MAKE) -C googletest all
	$(MAKE) -C GRAMCommon all
	$(MAKE) -C Earth all
endif

Mars:
ifneq ($(wildcard Mars),)
	$(MAKE) -C googletest all
	$(MAKE) -C GRAMCommon all
	$(MAKE) -C Mars all
endif

Jupiter:
ifneq ($(wildcard Jupiter),)
	$(MAKE) -C googletest all
	$(MAKE) -C GRAMCommon all
	$(MAKE) -C Jupiter all
endif

Neptune:
ifneq ($(wildcard Neptune),)
	$(MAKE) -C googletest all
	$(MAKE) -C GRAMCommon all
	$(MAKE) -C Neptune all
endif

Saturn:
ifneq ($(wildcard Saturn),)
	$(MAKE) -C googletest all
	$(MAKE) -C GRAMCommon all
	$(MAKE) -C Saturn all
endif

Titan:
ifneq ($(wildcard Titan),)
	$(MAKE) -C googletest all
	$(MAKE) -C GRAMCommon all
	$(MAKE) -C Titan all
endif

Uranus:
ifneq ($(wildcard Uranus),)
	$(MAKE) -C googletest all
	$(MAKE) -C GRAMCommon all
	$(MAKE) -C Uranus all
endif


clean:
	$(MAKE) -C googletest clean
	$(MAKE) -C GRAMCommon clean
ifneq ($(wildcard Venus),)
	$(MAKE) -C Venus clean
endif
ifneq ($(wildcard Earth),)
	$(MAKE) -C Earth clean
endif
ifneq ($(wildcard Mars),)
	$(MAKE) -C Mars clean
endif
ifneq ($(wildcard Jupiter),)
	$(MAKE) -C Jupiter clean
endif
ifneq ($(wildcard Neptune),)
	$(MAKE) -C Neptune clean
endif
ifneq ($(wildcard Saturn),)
	$(MAKE) -C Saturn clean
endif
ifneq ($(wildcard Titan),)
	$(MAKE) -C Titan clean
endif
ifneq ($(wildcard Uranus),)
	$(MAKE) -C Uranus clean
endif
	$(MAKE) -C GRAM clean
	
