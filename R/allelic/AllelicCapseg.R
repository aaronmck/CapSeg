
AllelicCapseg = function( capseg.probe.fn, capseg.seg.fn, germline.het.fn, SID, base.output.dir, min.seg.size, drop.x, drop.y,  verbose=FALSE )
{
   RESULTS.DIR=file.path(base.output.dir, "results" ); dir.create(RESULTS.DIR, recurs=TRUE, showWarnings=FALSE)
   plots.dir=file.path(base.output.dir, "plots", SID ); dir.create(plots.dir, recurs=TRUE, showWarnings=FALSE)

   result_FN = file.path(base.output.dir, paste(SID, ".AllelicCapseg.rds", sep="") )

 ## tmp.names is needed because the ReadGladMat parser requires the sample name corresponding to the rows it will extract so that it can pull single samples out of merged seg files.  However sometimes hyphens get mutated into dots and dumb stuff like that.  So tmp.name is what is actually in the seg file.
   tmp.name <- read.delim(capseg.seg.fn, stringsAsFactors=FALSE, check.names=FALSE)$Sample[1]
   seg.dat <- as.data.frame(ReadGladMat(capseg.seg.fn, sample.name=tmp.name, glad.log=TRUE, drop.x=drop.x, drop.y=drop.y, type="capseg", verbose=verbose)[[1]], stringsAsFactors=FALSE)

   if( !file.exists(result_FN))
   {
## Extract capture data: build list structures with germline het calls and capseg probe intensities according to the seg.dat created above and input here.
      cap.dat <- ExtractCaptureDat(capseg.probe.fn, seg.dat, germline.het.fn, drop.x=drop.x, drop.y=drop.y, verbose=verbose)

## Join small segments
#    merged.cap.res <- JoinSmallSegsExtreme(cap.dat, min.seg.size, verbose=verbose)

## Run model-fitting algorithm
      iams.res <- cap.dat
      iams.res[["em.fit"]] <- CaptureHscrSegFit(cap.dat[["as.res"]][["h.seg.dat"]], tol=1e-5, verbose=verbose)

      saveRDS( iams.res, file=result_FN)
   }
   else
   {
      iams.res = readRDS(result_FN)
   }

## Save output
   out.tab <- AbsolutePostProcess(iams.res, seg.dat)
   write.tab(out.tab, file=file.path(RESULTS.DIR, paste(SID, ".tsv", sep="")) )

## Plotting ##
   # images are automatically saved in RESULTS.DIR / plots / subdir
   PlotAllCapture(iams.res, save=TRUE, plots.dir )
}


##  TODO - remove dependency on these 'convenience functions'
ScarterPar = function()
{
   par(las=1)
   par(bty='n')
   par(cex=0.5)
}

trim = function(str)
{
    str = gsub("^\\s+", "", str);
    str = gsub("\\s+$", "", str);
    str
}

write.tab = function(..., sep = "\t", quote = F, row.names = F)
{
   write.table(..., sep = sep, quote = quote, row.names = row.names)
}



