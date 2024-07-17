#--------------------------------------------------------------------------
# Spherical reader

read_sph <- function(filename) {
   
  sphvoids <- read.table(filename, header=FALSE)
  colnames(sphvoids) <- c("id", "r", "x", "y", "z", "delta")

  nsph <- length(sphvoids$id)
  print(c("Number of spheres: ", nsph))

  return(sphvoids)

}

#--------------------------------------------------------------------------
# Popcorn reader

read_pop <- function(filename){

    con  <- file(filename, open="r")
    npop <- as.numeric(scan(con, what="", nlines=1, quiet=TRUE)) # an integer
    print(c("Number of popcorns: ", npop))

    for(i in 1:npop){

        head_pop <- as.numeric(scan(con, what="", nlines=1, quiet=TRUE)) # line that contains 5 numbers
        id    = head_pop[1] # an integer
        nmem  = head_pop[2] # an integer
        vol   = head_pop[3] # a float
        npart = head_pop[4] # an integer
        flag  = head_pop[5] # an integer

        x=NULL; y=NULL; z=NULL; r=NULL; fvol=NULL; level=NULL
        if(nmem>0){
            for(j in 1:nmem){
                popmem <- as.numeric(scan(con, what="", nlines=1, quiet=TRUE)) # line that contains 6 numbers
                x = c(x,popmem[1]) # adding a float
                y = c(y,popmem[2]) # adding a float
                z = c(z,popmem[3]) # adding a float
                r = c(r,popmem[4]) # adding a float
                fvol = c(fvol,popmem[5]) # adding a float
                level = c(level,popmem[6]) # adding an integer	
            }
        }

        if(npart>0){
            for(j in 1:npart){
                ttt <- as.numeric(scan(con, what="", nlines=1, quiet=TRUE))
		# integers that are not stored 
            }
        }
 
        pop = data.frame(x=x, y=y, z=z, r=r)
        df = list(id=id, nmem=nmem, vol=vol, npart=npart, pop=pop)

        if(i==1){
          popcorn=list(df)
        } else {
          popcorn[[i]]=df
        }

    }

    close(con)
    return(popcorn)

}

#--------------------------------------------------------------------------
# Effective popcorn spheres

pop2sph <- function(popcorn){

    npop <- length(popcorn)
    print(c("Number of pop-spheres: ", npop))
    
    id = NULL; x = NULL; y = NULL; z = NULL; r = NULL; vol = NULL
    for(i in 1:npop){
        id  <- c(id,  popcorn[[i]]$id)
        x   <- c(x,   popcorn[[i]]$pop$x[1])
        y   <- c(y,   popcorn[[i]]$pop$y[1])
        z   <- c(z,   popcorn[[i]]$pop$z[1])
        r   <- c(r,   popcorn[[i]]$pop$r[1])
        vol <- c(vol, popcorn[[i]]$vol)
        
        reff <- (vol*3/(4*pi))**(1/3)
    }
    
    effvoids <- data.frame(id, reff, x, y, z)
    colnames(effvoids) <- c("id", "r", "x", "y", "z")

    return(effvoids)

}
