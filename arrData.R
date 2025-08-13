dstORIG <- read.table("/home/dmorina/Insync/2102177@uab.cat/OneDrive Biz/Docència/UAB/2025-2026/Primer Semestre/Taller de Modelització/Data/DST_index.txt", sep=";", header=TRUE)

dstORIG$CDATE <- paste(dstORIG$DATE, dstORIG$TIME, sep=" ")
dstORIG$CDATE <- as.POSIXlt(dstORIG$CDATE, tz="UTC", format="%Y-%m-%d %H:%M:%OS")

dstORIG$DATE <- NULL
dstORIG$TIME <- NULL
dstORIG$DOY  <- NULL

dstORIG <- dstORIG[, c(2, 1)]

for (i in 1:40)
{
  dstORIG[paste0("STORM",i)] <- 0
  eval(parse(text=paste0("dstORIG$STORM", i, "[dstORIG$DST < ", -10*i, "] <- 1")))
  for (j in 2:dim(dstORIG)[1])
  {
    while(dstORIG[j, paste0("STORM", i)] == 1)
    {
      if (dstORIG$DST[j-1] < -10*i & dstORIG$DST[j] < -10*i)  dstORIG[j, paste0("STORM", i)] <- 0
      j <- j + 1
    }
  }
}

rm(i, j)
save.image("/home/dmorina/Insync/2102177@uab.cat/OneDrive Biz/Docència/UAB/2025-2026/Primer Semestre/Taller de Modelització/Data/dstORIG2025.RData")

#### Time bewtween two consecutive storms with the same threshold (including storms with less than 2 days of separation)
load("/home/dmorina/Insync/2102177@uab.cat/OneDrive Biz/Docència/UAB/2025-2026/Primer Semestre/Taller de Modelització/Data/dstORIG2025.RData")
for (i in 1:40)
{
  assign(paste0("diff", i), vector())
  assign(paste0("storm", i), dstORIG[dstORIG[,paste0("STORM", i)] == 1, ])
  for (j in 2:eval(parse(text=paste0("(dim(storm",i,")[1])"))))
  {
    eval(parse(text=paste0("diff",i,"[",j,"]",
                           "<-difftime(storm",i,"$CDATE[",j,"], storm",i,"$CDATE[",j-1,"], units='days')")))
  }
}
for (i in 1:40)
{
  eval(parse(text=paste0("length(diff",i,") <- length(diff1)")))
}

diff.times <- data.frame(diff1,diff2,diff3,diff4,diff5,diff6,
                         diff7,diff8,diff9,diff10,diff11,diff12,
                         diff13,diff14,diff15,diff16,diff17,diff18,
                         diff19,diff20,diff21,diff22,diff23,diff24,
                         diff25,diff26,diff27,diff28,diff29,diff30,
                         diff31,diff32,diff33,diff34,diff35,diff36,
                         diff37,diff38,diff39,diff40)

rm(list=setdiff(ls(), "diff.times"))
save(diff.times, file="/home/dmorina/Insync/2102177@uab.cat/OneDrive Biz/Docència/UAB/2025-2026/Primer Semestre/Taller de Modelització/Data/diff.timesORIG2025.RData")

### Remove storms with less than 2 days of separation
load("/home/dmorina/Insync/2102177@uab.cat/OneDrive Biz/Docència/UAB/2025-2026/Primer Semestre/Taller de Modelització/Data/dstORIG2025.RData")
dst <- dstORIG
pr  <- dst
pr2 <- 0
for (k in 1:40)
{
  pr  <- dst
  pr2 <- 0
  while(min(pr2)<2)
  {
    for (i in 1:dim(dst)[1])
    {
      if (eval(parse(text=paste0("dst$STORM", k, "[i]==1"))))
      {
        j <- 1
        while(eval(parse(text=paste0("dst$STORM", k, "[i+j]==0 & (i+j)<=dim(dst)[1]"))))
        {   
          j <- j + 1
        }
        if ((i+j) <= dim(dst)[1])
        {
          if (difftime(dst$CDATE[i+j], dst$CDATE[i], units="days")<2) 
          {
            eval(parse(text=paste0("dst$STORM", k, "[i+j] <- 0")))
          }
        }
      }
    }
    pr  <- dst[eval(parse(text=paste0("dst$STORM", k, "==1"))), ]
    pr2 <- diff(pr$CDATE)
    units(pr2) <- "days"
  }
}

rm(list=setdiff(ls(), "dst"))
save.image("/home/dmorina/Insync/2102177@uab.cat/OneDrive Biz/Docència/UAB/2025-2026/Primer Semestre/Taller de Modelització/Data/dstDEF2025.RData")

#### Time bewtween two consecutive storms with the same threshold (excluding storms with less than 2 days of separation)
load("/home/dmorina/Insync/2102177@uab.cat/OneDrive Biz/Docència/UAB/2025-2026/Primer Semestre/Taller de Modelització/Data/dstDEF2025.RData")
for (i in 1:40)
{
  assign(paste0("diff", i), vector())
  assign(paste0("storm", i), dst[dst[,paste0("STORM", i)] == 1, ])
  for (j in 2:eval(parse(text=paste0("(dim(storm",i,")[1])"))))
  {
    eval(parse(text=paste0("diff",i,"[",j,"]",
                           "<-difftime(storm",i,"$CDATE[",j,"], storm",i,"$CDATE[",j-1,"], units='days')")))
  }
}
for (i in 1:40)
{
  eval(parse(text=paste0("length(diff",i,") <- length(diff1)")))
}

diff.times <- data.frame(diff1,diff2,diff3,diff4,diff5,diff6,
                         diff7,diff8,diff9,diff10,diff11,diff12,
                         diff13,diff14,diff15,diff16,diff17,diff18,
                         diff19,diff20,diff21,diff22,diff23,diff24,
                         diff25,diff26,diff27,diff28,diff29,diff30,
                         diff31,diff32,diff33,diff34,diff35,diff36,
                         diff37,diff38,diff39,diff40)

rm(list=setdiff(ls(), "diff.times"))
save(diff.times, file="/home/dmorina/Insync/2102177@uab.cat/OneDrive Biz/Docència/UAB/2025-2026/Primer Semestre/Taller de Modelització/Data/diff.timesDEF2025.RData")