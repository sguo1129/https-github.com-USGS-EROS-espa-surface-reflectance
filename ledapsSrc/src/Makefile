
MODLIST = lndpm lndcal lndsr lndsrbm 6sV-1.0B

all:
	@for i in $(MODLIST); do \
	(cd $$i; make; cd ..) \
	done

install:
	@for i in $(MODLIST); do \
	(cd $$i; make install; cd ..) \
	done
	install -m 755 ../bin/do_ledaps.py $(BIN)

clean:
	@for i in $(MODLIST); do \
	(cd $$i; make clean; cd ..) \
	done
