  /*  printf("%s %s totint %g totint_bulk %g totprimint %g totprimintbulk %g"
                     "totsecint %g totcomptonint %g totcomptonintbulk %g\n"
               "totraylfotint %g totraylfotint_bulk %g \n"
               "totfotraylint %g totfotraylint_bulk %g  \n",
            el, lin, totint, totint_bulk, totprimint, totprimint_bulk,
                         totsecint, totcomptonint, totcomptonint_bulk,
                  totraylfotint, totraylfotint_bulk, totfotraylint, totfotraylint_bulk);
	printf("rxi old is %g rxi scatter is %g \n", 
		(totprimint+totsecint)/(totint_bulk-totcomptonint_bulk-totraylfotint_bulk-
		   totfotraylint_bulk), totint/totint_bulk); /* */
    rxi_calc[j]=totint/totint_bulk;
