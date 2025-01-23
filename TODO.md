# TODO list

* 202501: Check number of conditions, if > 8, change the default
  palette from 'Dark2' to something with >=x.
* 202406: ontology stuff: I want to redo the clusterProfiler run/write
  to make explicit the distinction between the over representation
  analyses (high/low/up/down) vs gene set enrichment (distribution of
  values (expression/fc/etc)).  I started this process in
  write_cp_data but it needs fleshing out.
* 202309: write_expt(), check the design rank before calling
  simple_varpart() and/or add some logic to simple_varpart() so that
  it a) does not trip over non-filtered data and b) nor non-qr passing
  models.
* 202309: Finish creating methods and removing all the various
  if(class()) stanzas.

