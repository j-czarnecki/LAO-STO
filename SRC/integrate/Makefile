.PHONY: test clean profile

test:
	$(MAKE) -C src all
	$(MAKE) -C test all
	$(MAKE) -C test run
profile:
	$(MAKE) -f test/profile.mk profile
	$(MAKE) -f test/profile.mk run
clean:
	$(MAKE) -C src clean
	$(MAKE) -C test clean
	$(MAKE) -f test/profile.mk clean