## Coevolution with a seed bank - analysis of sporulation in evolved lines
library(renv)
# # # renv::init()
# renv::restore()
library(here)
#analysis of flow-cytometry population data using Karava's analysis.R as reference.
source(here("src/FCM_functions.R"))
set.seed(1)
#point to data
   all.folders <- 
   list.dirs(here("data/FCM/"),recursive = F)

for (current_folder in all.folders){
   
   # day is an experimental batch
   day <- str_extract(current_folder, "FCM/.*$") %>% 
      str_remove("FCM/")

   folders <- list.dirs(current_folder)
   # remove the folder of folders
   if(length(folders) > 1){
      folders <- folders[!str_detect(folders, paste0(day,"$"))]
   }
   # place "xIPTG" first" to be used in making of stran model
   tmp <- c(folders[grep("xIPTG", folders)], folders[grep("xIPTG", folders, invert = TRUE)])
   folders <- tmp
   rm(tmp)
   
   # Make directories to store the data
   if (! dir.exists(here("fig/gate_plots", day))){
      dir.create(here("fig/gate_plots", day), recursive = T)
   }
   
   if (! dir.exists(here("data/output", day))){
      dir.create(here("data/output", day), recursive = T)
   }
   #fields to parse from the .fcs file name (separated by underscore and/or hyphen)
   sample.var <- c("strain","media","treat","dilution","well","colony" ,"rep","xt","num") 
   
   # the following loop goes over each replicated samples and does the actual anaysis
   for (folder in folders){
      #### Load data, sample set ####
      fcsset <- flowCreateFlowSet(filepath = folder, sample_variables = sample.var,
                                  transformation = FALSE,separators = "[-_\\.]")
      #transform with arcsine, recomendded by Karava et al.
      fcsset <- Transform.Novocyte(fcsset)
      
      #data frame to collect stats
      df.stats <- fcsset%>%
         flowFcsToDf(.)%>%
         select(.,sample.var)%>%
         distinct()
      df.stats$exp <- day
      df.stats$total.events <- NA
      df.stats$volume.nL <- NA
      df.stats$limit.events <- NA
      df.stats$limit.volume.uL <- NA
      df.stats$singlets <- NA
      df.stats$neg.rmv <- NA
      df.stats$noise.cutoff <- NA
      
      
      #### Gating for singlets with flowStats ####
      
      # The gate function needs to be applied to each sample separately
      # get number of samples
      n.sample <- nrow(fcsset@phenoData@data)
      
      
      
      #initialise list to store singlet plots
      plot.list <- vector('list', n.sample)
      
      for (i in 1:n.sample){
         # collect metadata for stats table
         df.stats$volume.nL[i] <- as.numeric(fcsset[[i]]@description$`$VOL`)
         df.stats$total.events[i] <- as.numeric(fcsset[[i]]@description$`$TOT`)
         df.stats$limit.volume.uL[i] <- as.numeric(fcsset[[i]]@description$`#NCVolumeLimits`)
         df.stats$limit.events[i] <- as.numeric(fcsset[[i]]@description$`#NCEventsLimits`)
         
         
         singlet_gate <- gate_singlet(fcsset[[i]], area = "asinh.FSC.A", height = "asinh.FSC.H", filterId = "Singlets",wider_gate = TRUE )
         df.stats$singlets[i] <- summary(filter(fcsset[[i]],singlet_gate))@true
         
         #plot gate
         id <- fcsset[[i]]@description$GUID
         plot.list[[i]] <-
            as.ggplot(
               ggcyto(fcsset[[i]], aes(x = `asinh.FSC.A`, y =  `asinh.FSC.H`))+
                  geom_hex(bins = 500)+
                  geom_gate(singlet_gate, size=0.01)+
                  geom_stats(adjust=c(0.2,0.75))+
                  theme_cowplot(font_size = 8)+
                  scale_y_continuous(labels =  function(x) format(x, scientific = TRUE))+
                  scale_x_continuous(labels =  function(x) format(x, scientific = TRUE))+
                  facet_null()+
                  ggtitle(df.stats$well[i])
               
            )
         
         #apply gate
         fcsset[[i]] <- fcsset[[i]] %>%
            Subset(singlet_gate)
         
         
         # filter negatives
         neg.gate <- rectangleGate("asinh.BL1.A" = c(0, 15), "asinh.FSC.A" = c(0, 15),"asinh.SSC.A" = c(0, 15))
         df.stats$neg.rmv[i] <-  df.stats$singlets[i]- summary(filter(fcsset[[i]],neg.gate))@true
         #apply gate
         fcsset[[i]] <- fcsset[[i]] %>%
            Subset(neg.gate) # remove negative
         
         #find cutoff to remove noise by FSC
         noise<- rangeGate(x = fcsset[[i]], "asinh.FSC.A", alpha=0.5,sd=3, absolute = F)
         df.stats$noise.cutoff[i] <- noise@min[[1]]
         
      }
      
      #save plot
      ggsave2(filename = paste0(fcsset[[i]]@description$`$SRC`,".pdf" ), 
              plot = plot_grid(plotlist = plot.list),
              path = here("fig/gate_plots", day))
      
      
      # plotting noise gate
      
      noise.plot <-   
         ggcyto(fcsset,aes(asinh.FSC.A))+
         geom_density(fill="grey80")+
         geom_gate(noise)+
         
         facet_wrap(~well)+
         theme_cowplot()+
         panel_border(color = "black", size=.2)+
         geom_stats(x=Inf, y=Inf, hjust=1.1, vjust=1.1)+
         xlab("asinh.FSC.A")
      
      ggsave2(filename = paste0("noise_",fcsset[[i]]@description$`$SRC`,".pdf" ), 
              plot = noise.plot,
              path = here("fig/gate_plots/", day))
      
      #### transform to dataframe ####
      df.set <- fcsset%>%
         flowFcsToDf(.)
      # no. of events before noise filter
      df.stats <- 
         df.set %>%
         group_by(well) %>%
         summarise(noisy.events=n())%>%
         full_join(df.stats,.)
      
      
      
      # remove noise
      clean.df <- df.set[0,]
      wells <- levels(as.factor(df.stats$well))
      for (wl in wells){
         clean.df <- rbind(clean.df,
                           df.set %>%
                              dplyr::filter(well==wl)%>%
                              dplyr::filter(asinh.FSC.A>df.stats$noise.cutoff[df.stats$well==wl]))
      }
      
      
      df.stats <- 
         clean.df %>%
         group_by(well) %>%
         summarise(clean.events=n())%>%
         full_join(df.stats,.)
      
      complot <- 
         ggcyto(fcsset,aes(asinh.FSC.A, asinh.BL1.A))+
         geom_hex(bins=300)+
         geom_gate(noise)+
         geom_stats(x=Inf, y=Inf, hjust=1, vjust=1)+
         # geom_density2d(col = "black",  size = 0.1, alpha=0.5) +
         facet_wrap(~well)+
         theme_cowplot()+
         panel_border()+
         ylab("asinh.SYBR-green.A")+
         xlab("asinh.FSC.A")
      
      ggsave2(filename = paste0("scatterNoise_",fcsset[[i]]@description$`$SRC`,".png" ), 
              plot = complot,
              path = here("fig/gate_plots", day))
      #### predict centers of sub-populations ####
      
      # I pre-compiled a model for Bacillus subtilis TS01
      # base::load("data/cluster_models/TS01_cluster_model.Rdata")
      # # use this to generate prediction model from data frame of events
      # df.mix <- clean.df%>%
      #       dplyr::select(asinh.FSC.A,asinh.BL1.A)%>%
      #       as.matrix() %>%
      #       mclust::Mclust(data = ., G = 2) #I  have only 2 clusters
      # save(df.mix, file="data/cluster_models/name_cluster_model.Rdata")
      
      current.strain <- paste(df.set$strain[1],day, sep = "_")
      models <- list.files("data/cluster_models/", full.names = T)
      
      
      if (!any(grepl(current.strain,models))&
          !grepl("blank",current.strain)){
         # generate cluster prediction model for strain from current set
         df.mix <- clean.df%>%
            # dplyr::filter(phage=="noPHI")%>%
            dplyr::select(asinh.FSC.A,asinh.BL1.A)%>%
            as.matrix() %>%
            mclust::Mclust(data = ., G = 2) #I  have only 2 clusters
         save(df.mix, file=paste0("data/cluster_models/",current.strain,"_cluster_model.Rdata")) 
         model.use <- paste("new for ", current.strain)
      } else{
         if (grepl("blank",current.strain)){
            base::load("data/cluster_models/TS01_cluster_model.Rdata")
            model.use <- "TS01"
         } else{
            # load previously compiled prediction model for strain
            model.use <- models[grep(current.strain,models)]
            base::load(model.use)
         }
      }
      df.stats$clust.model <- model.use
      
      # getting centers for visualization and export
      centers.list.df <- t(df.mix$parameters$mean)
      
      # assigning cluster to population based on SYBR fluoresence (BL1)
      # higher SYBR => veg cell
      pop.tbl <- data.frame(cluster=c(1,2), pop=NA)
      pop.tbl$pop[which.max(centers.list.df[,"asinh.BL1.A"])] <- "veg"
      pop.tbl$pop[which.min(centers.list.df[,"asinh.BL1.A"])] <- "spore"
      
      ### Prediction of clusters for all samples ####
      cluster.predict <- clean.df %>%
         dplyr::select( asinh.FSC.A, asinh.BL1.A) %>%
         predict.Mclust(df.mix,.)
      
      clean.df$pop <- sapply(cluster.predict$classification,  function(x) {ifelse(x==pop.tbl$cluster[1], pop.tbl$pop[1],pop.tbl$pop[2])})
      
      
      
      #### Get the quantities ####
      
      clust.count <- data.frame(clean.df, cluster = cluster.predict$classification) %>%
         group_by(well, cluster) %>%
         summarize(count = n()) 
      for(i in pop.tbl$cluster){
         clust.count$cluster[clust.count$cluster %in%  pop.tbl$cluster[i]] <-  pop.tbl$pop[i]
      }
      
      df.stats <-    
         spread(clust.count ,cluster, count)%>%
         full_join(df.stats,.)
      

      # calculate concentrations based on volume and dilution
      df.stats$spore.ml <- as.numeric(sapply(strsplit(df.stats$dilution,"x"), "[[", 2))* #dilution factor
         df.stats$spore/(df.stats$volume.nL/1e6)
      df.stats$veg.ml <- as.numeric(sapply(strsplit(df.stats$dilution,"x"), "[[", 2))* #dilution factor
         df.stats$veg/(df.stats$volume.nL/1e6)  
      
      # write results to file
      write_csv(df.stats,
                path = here("data/output/",day,paste0(fcsset[[i]]@description$`$SRC`,".csv")))
      
      
      
      # plot clusters
      p <-
         ggplot(clean.df, aes(asinh.FSC.A, asinh.BL1.A)) +
         geom_point(aes(color = pop), size=.1)+
         # geom_hex(aes(fill = pop), bins = 300) + # ,alpha=..ncount.. #order= ?
         geom_density2d(color="black",  size = 0.1) +
         geom_text(data = df.stats, aes(label = paste0(veg+spore, " events")), 
                   x=8, y = 3, hjust = 0, vjust = 0,
                   size=4)+
         geom_text(data = df.stats, aes(label = paste0(round(100*veg/(veg+spore),1), "% veg")), 
                   x=8, y = 1.5, hjust = 0, vjust = 0,
                   color = scales::hue_pal()(2)[2], size=4)+
         geom_text(data = df.stats, aes(label = paste0(round(100*spore/(veg+spore),1), "% spore")), 
                   x=8, y = 0, hjust = 0, vjust = 0,
                   color = scales::hue_pal()(2)[1], size=4)+
         theme_bw()+ 
         geom_point(aes(centers.list.df[1, 1], centers.list.df[1, 2]), col = "blue", size = 1) +
         geom_point(aes(centers.list.df[2, 1], centers.list.df[2, 2]), col = "blue", size = 1) +
         facet_wrap(~well)+
         ylab("asinh.SYBR-green.A")+
         ylim(0,15)+
         xlim(8,15)+
         guides(color = guide_legend(override.aes = list(size=3)))
      
      ggsave2(filename = paste0("cluster_",fcsset[[i]]@description$`$SRC`,".png" ), 
              plot = p,
              path = here("fig/gate_plots", day))
      
      print(paste("done",folder))
   } #folder loop
   
} # experimental batch loop     






