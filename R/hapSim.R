
###########################################
# generate genotype data by gene dropping #
###########################################

#######################
# recode the pedigree #
#######################

pedRecode <- function(ped,ids){
# ped: pedigree (id, sire, dam,...) or  (id, generation, sire, dam,...)
#    missing values in sire or dam represented by 0 or NA
# ids: IDs of interest if not missing
# output: id=1,2,...
   ped<- as.data.frame(ped)
      pedSave<- ped
   if(is.null(ped$id)){
      stop("'id' missing...")
   }else ped$id<- trim(ped$id)
   if(is.null(ped$sire)){
      stop("'sire' missing...")
   }else ped$sire<- trim(ped$sire)
   if(is.null(ped$dam)){
      stop("'dam' missing...")
   }else ped$dam<- trim(ped$dam)

   idx<- is.na(ped$id) | ped$id==0 | ped$id=="0"
   if(any(idx)){
      print(ped[idx,])
      cat("   Above individuals with N/A IDs were removed.\n")
      ped<- ped[!idx,]
   }

   if(!missing(ids)){# discard irrelevant ones
      ids<- trim(ids)
      if(length(ids)==0) stop("IDs not correctly specified.")
      idx<- !is.element(ids,ped$id)
      if(any(idx)){
         print(ids[idx])
         stop("Check the above IDs; out of range.")
      }

      idTmp<- ids
      idx<- rep(FALSE,nrow(ped))
      while(1){
         idxTmp<- is.element(ped$id,idTmp)
         idx<- idx | idxTmp
         idx<- idx | is.element(ped$id,ped$sire[idxTmp]) |
                     is.element(ped$id,ped$dam[idxTmp])
         idTmp<- ped$id[idx]

         if(sum(idxTmp)==sum(idx)) break
      }
      ped<- ped[idx,]
      rm(idTmp,idx,idxTmp)
   }

   if(is.null(ped$generation))
      ped<- pedRecode.0(ped)
   ids<- paste(ped$generation,ped$id,sep="~")
   uids<- unique(ids)
      idx<- match(uids,ids)
   if(length(uids)<length(ids)){
      cat("   The following samples are repeated:\a\n")
      print(ped[-idx,c("id","generation","sire","dam")])
      cat("   Repeated IDs were excluded!\n")
   }
   ped<- ped[idx,]
   rm(idx,ids)

   # new code
   ped$generation<- reorder(factor(ped$generation))
   ped$id<- reorder(factor(ped$id))
   ped$sire<- reorder(factor(ped$sire))
   ped$dam<- reorder(factor(ped$dam))
   ped<- ped[order(ped$generation,ped$id,ped$sire,ped$dam),]
   idd<- data.frame(index=c(0,0),id=c(NA,0))
      idd<- rbind(idd,data.frame(index=1:nrow(ped),id=ped$id))
   ### recode here

   # recode IDs
   if(is.null(ped$sex)){
      idx<- match(ped$sire,idd$id)
         idx[is.na(idx)]<- 1
      sire<- idd$index[idx]
      idx<- match(ped$dam,idd$id)
         idx[is.na(idx)]<- 1
      dam<- idd$index[idx]
      ped<- data.frame(id=idd$index[-c(1:2)],
                       sire=sire,
                       dam=dam,
                       generation=ped$generation,
                       old.id=ped$id)
   }else{
      idx<- match(ped$sire,idd$id)
         idx[is.na(idx)]<- 1
      sire<- idd$index[idx]
      idx<- match(ped$dam,idd$id)
         idx[is.na(idx)]<- 1
      dam<- idd$index[idx]
      ped<- data.frame(id=idd$index[-c(1:2)],
                       sire=sire,
                       dam=dam,
                       sex=ped$sex,
                       generation=ped$generation,
                       old.id=ped$id)
      ii<- match(ped$sire,ped$id)
         ii<- ii[!is.na(ii)]
      if(length(ii)>0) if(any(ped$sex[ii]!=1 & ped$sex[ii]!="M")){
       #  cat("   Suppose 1 or \"M\" stands for male...\n")
       #  print(ped[ii,][ped$sex[ii]!=1 & ped$sex[ii]!="M",])
       #  stop("   Check above sires for errors in sex.")
      }
      jj<- match(ped$dam,ped$id)
         jj<- jj[!is.na(jj)]
      if(length(jj)>0) if(any(ped$sex[jj]==1 || ped$sex[jj]=="M")){
         # cat("   Suppose !0 or \"!M\" stands for female...\n")
         # print(ped[jj,][ped$sex[jj]==1 || ped$sex[jj]=="M",])
         # stop("   Check above dams for errors in sex.")
       }
   }
   idx<- (ped$sire > ped$id) | (ped$dam > ped$id)
   if(any(idx)){
      idx<- match(ped$old[idx],pedSave$id)
      print(pedSave[idx,][1:min(sum(idx),3),])
      cat("... ...\n")
      stop("Pedigree might not correct. Check the above for errors.")
   }
   rownames(ped)<- 1:nrow(ped)

   ped
}

# create "generation"
pedRecode.0<- function(ped){
# ped: data frame (id,sire,dam,...)
   ped<- as.data.frame(ped)
      ped$generation<- NULL
   if(is.null(ped$id)){
      stop("id missing...")
   }else ped$id<- trim(ped$id)
   if(is.null(ped$sire)){
      stop("sire missing...")
   }else ped$sire<- trim(ped$sire)
   if(is.null(ped$dam)){
      stop("dam missing...")
   }else ped$dam<- trim(ped$dam)

   idx<- is.na(ped$id) | ped$id==0 | ped$id=="0"
   if(any(idx)){
      print(ped[idx,])
      cat("   above individuals with N/A IDs were removed.\n")
      ped<- ped[!idx,]
   }

   nr<- nrow(ped)
   ii<- matrix(TRUE,nrow=1,ncol=nr)
   ii0<- ii
   while(1){
      ii0<- is.element(ped$id,ped$sire[ii0]) |
            is.element(ped$id,ped$dam[ii0])
      if(any(ii0)){
         ii<- rbind(ii0,ii)
      }else break
   }
# no offspring or parents
#   idx<- is.element(ped$id,ped$sire) |
#         is.element(ped$id,ped$dam)  |
#         is.element(ped$sire,ped$id) |
#         is.element(ped$dam,ped$id)
#   ii[1,]<- ii[1,] | !idx

   idx<- ii[1,]
   idx0<- idx
   jj<- 0
   out<- cbind(generation=jj,ped[idx,])
   if(nrow(ii)>1){
      for(n in 2:nrow(ii)){
         idx0<- idx0 | ii[n-1,]
         idx<- ii[n,] & !idx0
         if(any(idx)){
            jj<- jj+1
            out<- rbind(out,cbind(generation=jj,ped[idx,]))
         }
      }
   }
   jj<- match("id",colnames(out))
   out<- cbind(id=out[,jj],generation=out$generation,out[-c(1,jj)])
   rownames(out)<- 1:nrow(out)

   out$generation<- reorder(factor(out$generation))
   out$id<- reorder(factor(out$id))
   out$sire<- reorder(factor(out$sire))
   out$dam<- reorder(factor(out$dam))

   out<- out[order(out$generation,out$id,out$sire,out$dam),]
   out
}

#######################################
# generate haplotype or genotype data #
#######################################
# the follow code generate genotype data chrosomome by chromosome

.hapSim<- function(ped,gmap,ids,hap,method=c("Haldane","Kosambi"),recode.pedigree=FALSE,genotype=FALSE){
   if(recode.pedigree){
      idTmp<- ped$id
      ped<- pedRecode(ped)
      if(!missing(hap)){
         if(any(hap<0)) stop("hap: Use non-negative integers for alleles.")
         if(any(hap>2^14)) stop("hap: Too large integers.")
         idx<- ped$generation=="F0" | ped$generation=="0"
         if(!any(idx)) stop("check pedigree for errors in founders' generation.")
            idx<- match(ped$old.id[idx],idTmp[1:nrow(hap)])
         if(length(idx)>nrow(hap)) stop("haplotypes are not specified for all founders?")
         hap<- hap[idx,]
      }
   }
   if(!is.data.frame(gmap) || any(!is.element(c("snp","chr","dist"),colnames(gmap)))){
      stop("genetic map should be a data frame (snp,chr,dist,...).")
   }
   if(!missing(ids)){
      ids<- trim(ids)
      if(length(ids)==0) stop("ids not correctly specified.")
      ii<- match(ids,ped$old.id)
      if(any(is.na(ii))) stop("check ids for error.")
   }else ii<- 1:length(ped$old.id)

   gmap$chr<- reorder(factor(gmap$chr))
      ord<- order(gmap$chr,gmap$dist)
      gmap<- gmap[ord,]
   if(is.null(gmap$recRate)){
      if(!is.numeric(gmap$dist)){
         stop("gmap$dist should be numeric.")
      }
      method<- match.arg(method)
      dstTmp<- diff(gmap$dist)
         dstTmp[dstTmp<0]<- 0
         dstTmp<- c(0,dstTmp)/100
      gmap$recRate<- mappingFuncInv(dstTmp,method=method)
   }
   chr<- unique(gmap$chr)
      chr<- reorder(chr)
      chr<- sort(chr)
   nchr<- length(chr)
   out<- NULL
   for(n in 1:nchr){
      seed<- runif(1,min=0,max=2^31-1)
         seed<- round(seed,0)
      tmp<- .hapSim0(ped,gmap,chr[n],hap,seed,genotype)
      out<- cbind(out,tmp)
   }
   out<- out[ii,];
      out<- as.matrix(out)
   rownames(out)<- ped$old.id[ii]
   idx<- order(ord)
   if(genotype){
      as.matrix(out[,idx])
   }else{
      idx<- rbind(idx*2-1,idx*2)
      as.matrix(out[,idx])
   }
}

.hapSim0<- function(pedd,gmap,chr,hap,seed,genotype){
   pedd<- pedd[,c("id","sire","dam","sex")]
   if(is.numeric(pedd[,"sex"])){
      pedd[,"sex"]<- pedd[,"sex"]==1
   }else{
      pedd[,"sex"]<- pedd[,"sex"]=="M" | pedd[,"sex"]=="Male" |
                     pedd[,"sex"]=="m" | pedd[,"sex"]=="male"
   }
   idx<- gmap$chr==chr
   rr<- gmap$recRate[idx]
   nr<- nrow(pedd)
   nc<- sum(idx)
   xchr<- FALSE
      if(chr=="x" || chr=="X") xchr<- TRUE
   if(missing(seed)) seed<- 0
   gdat<- matrix(-99,nrow=nr,ncol=2*nc)
   if(!missing(hap)){
      ninit<- nrow(hap)
      gdat[1:ninit,]<- hap[,rep(idx,rep(2,length(idx)))]
   }else{
      ninit<- 2
      if(xchr){
         gdat[1,]<- rep(c(0,1),nc)
      }else gdat[1,]<- rep(c(1,1),nc)
      gdat[2,]<- rep(c(2,2),nc)
   }

   out<- .C("rgdata2",
            gdata = as.integer(t(gdat)),
            nr = as.integer(nr),
            nc = as.integer(nc),
            ninit = as.integer(ninit),
            pedigree = as.integer(t(pedd)),
            recomb = as.double(rr),
            xchr = as.logical(xchr),
            seed = as.integer(seed),
            DUP = FALSE)$gdata
   out[out==-99]<- NA
   out<- matrix(out,nrow=nr,byrow=TRUE)
      storage.mode(out)<- "integer"

   if(genotype){
      oo<- out[,2*(1:nc)-1] + out[,2*(1:nc)]; oo<- as.matrix(oo)
      if(xchr){
         ii<- as.logical(pedd[,"sex"])
         if(any(out[ii,2*(1:nc)-1]!=0))
            stop("paternal wrong! check pedigree for sex errors.")
         if(any(out[!ii,2*(1:nc)-1]==0))
            stop("maternal wrong! check pedigree for sex errors.")
         oo[ii,]<- 2*oo[ii,]
      }
      oo<- oo-1
      storage.mode(oo)<- "integer"
      colnames(oo)<- gmap$snp[gmap$chr==chr]
      oo
   }else{
      out
   }
}

hapSim<- function(ped,gmap,ids,hap,method=c("Haldane","Kosambi"),recode.pedigree=FALSE){
# ped: recoded pedigree
# gmap: genetic map (snp,chr,recom,dist,...)
# hap: founders' haplotypes
# ids: only output data for individuals with ID ids
   if(missing(gmap)){
      gmap<- data.frame(snp="N",chr="N",dist=0)
      if(!missing(hap)){
         hap<- as.matrix(hap)
         if(ncol(hap) < 2) stop("hap should have at least 2 columns.")
         hap<- hap[,1:2]
      }
   }
   out<- .hapSim(ped = ped,
                gmap = gmap,
                ids = ids,
                hap = hap,
                method = method,
                recode.pedigree = recode.pedigree,
                genotype = FALSE)
   out
}

genoSim<- function(ped,gmap,ids,hap,method=c("Haldane","Kosambi"),recode.pedigree=FALSE){
# ped: recoded pedigree
# gmap: genetic map (snp,chr,recom,dist,...)
# output: individuals by SNPs
#   for each SNP 1--AA, 2--AB,3--BB
# ids: only output data for individuals with ID ids
   if(missing(gmap)){
      gmap<- data.frame(snp="N",chr="N",dist=0)
      if(!missing(hap)){
         hap<- as.matrix(hap)
         if(ncol(hap) < 2) stop("hap should have at least 2 columns.")
         hap<- hap[,1:2]
      }
   }
   out<- .hapSim(ped = ped,
                gmap = gmap,
                ids = ids,
                hap = hap,
                method = method,
                recode.pedigree = recode.pedigree,
                genotype = TRUE)
   out
}

