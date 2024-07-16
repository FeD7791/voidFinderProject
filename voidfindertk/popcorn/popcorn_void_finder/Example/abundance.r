source("read_voids.r")

# Reading the spherical voids catalogue:
sphvoids <- read_sph("sphvoids.dat")

# Reading the popcorn voids catalogue:
popvoids <- read_pop("popvoids.dat")

# Getting the effective popcorn spheres:
effvoids <- pop2sph(popvoids)

# Calculating radii histograms:
hs <- hist(sphvoids$r, plot=FALSE)
hp <- hist(effvoids$r, plot=FALSE)

# Plotting radii distribution:
plot( hs$mids, hs$counts
    , log="xy"
    , type="l", lwd=2, col="blue"
    , xlab=expression(R[v]~"[Mpc/h]")
    , ylab="Void counts" 
    )

lines( hp$mids, hp$counts
     , lwd=2, col="red"
     )

legend( "topright"
      , c("spherical", "popcorn")
      , lty=1, lwd=2, col=c("blue", "red")
      , bty="n" 
      )
