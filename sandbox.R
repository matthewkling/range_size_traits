library(dplyr)
library(tidyr)
library(ggplot2)

############ TRAITS ############

trr <- readRDS("E:/BIEN/traits_ranges/data/traits.rds")
tr <- trr %>%
      distinct() %>%
      select(family, genus, taxon, trait_name, trait_value) %>%
      na.omit() %>%
      group_by(family, genus, taxon, trait_name) %>%
      summarize(trait_value=mean(as.numeric(trait_value), na.rm=T)) %>%
      mutate(trait_name = gsub(" ", "_", trait_name)) %>%
      filter(trait_value != 0) %>%
      spread(trait_name, trait_value) %>%
      mutate(Leaf_NP_ratio = Leaf_Nmass / Leaf_Pmass) %>%
      gather(trait_name, trait_value, -taxon, -genus, -family) %>%
      na.omit()

# records per trait
ggplot(tr, aes(trait_name)) + geom_histogram()
frq <- as.data.frame(table(tr$trait_name))
tr <- filter(tr, trait_name %in% frq$Var1[frq$Freq>1000] &
                   trait_name != "Flowering_month") # only include widely recorded traits

# add habit data
h <- read.table("E:/BIEN/traits_ranges/data/Habit_Final.txt", sep="\t", header=T) %>%
      mutate(taxon=Accepted_name, habit=BIENHABIT) %>%
      select(taxon, habit)
tr <- left_join(tr, h)
frq <- as.data.frame(table(tr$habit))
tr <- filter(tr, habit %in% frq$Var1[frq$Freq>1000]) # only include widely recorded traits


# trait distributions by habit
p <- ggplot(tr, aes(habit, trait_value, color=habit, fill=habit)) +
      geom_boxplot(alpha=.5) +
      scale_y_log10() +
      facet_wrap(~trait_name, scales="free", ncol=4) +
      labs(title="Trait distributions by habit", y="trait value (log scale)") +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            legend.position="top")
ggsave(paste0("E:/BIEN/traits_ranges/v5/trait_by_habit_boxplot.png"), p, width=12, height=9, units="in")


# trait covariance
#tcd <- tr %>%
#      filter(is.finite(trait_value)) %>%
#      spread(trait_name, trait_value) %>%
#      dplyr::select(-taxon, -genus, -family) %>%
#      as.matrix() %>%
#      log10()
#pairs(tcd)
#corrplot::corrplot(cor(tcd, use="pairwise.complete.obs"))
#pc <- prcomp(na.omit(tcd)) # NO species has data on all 7 of the most common traits




############ range/rairity ############

load("E:/BIEN/traits_ranges/data/allNewAreas.rdata")

rr <- Areas %>%
      mutate(taxon=Latin) %>%
      select(-Latin)

rd <- gather(rr, variable, value, -taxon)
ggplot(rd, aes(value)) +
      geom_density() +
      facet_wrap(~variable, scales="free") +
      scale_x_log10()

# range variable covariance (not informative due to structure of missing data)
#pairs(rr[,1:5])
#rrd <- log10(as.matrix(rr[,1:5]))
#pcrr <- prcomp(na.omit(rrd[,c(1,2,5)]))
#pcrr <- data.frame(taxon=rr$taxon, range=pcrr$x[,1])

rr <- select(rr, taxon, sampleSize, clippedHullArea, updArea)


########### join datasets ##########
d <- left_join(tr, rr) %>%
      gather(range_name, range_value, -family, -genus, -taxon, -trait_name, -trait_value, -habit) %>%
      #mutate(trait_name_break = gsub(" ", "\n", trait_name)) %>%
      na.omit()

#d <- filter(d, range_name %in% c("clippedHullArea", "sampleSize"))


######### basic scatterplots at 3 taxanomic levels #########
for(level in c("species", "genus", "family")){
      
      dl <- d
      
      title <- paste0(level, "-level")
      
      if(level %in% c("genus", "family")){
            dl <- dl %>%
                  group_by_(level, "trait_name", "range_name") %>%
                  summarize(trait_value=mean(trait_value),
                            range_value=mean(range_value))
      }
      
      p <- ggplot(na.omit(dl[dl$trait_value>=0 & dl$range_value>=0,]),
                  aes(trait_value, range_value)) +
            geom_point(size=1.5, alpha=.3) +
            geom_smooth(color="dodgerblue", se=F, size=1) +
            geom_smooth(color="red", se=F, size=1, method=lm) +
            facet_grid(range_name~trait_name, scales="free") +
            labs(title=paste0(title, ", raw\n"))
      #ggsave(paste0("E:/BIEN/traits_ranges/v5/scatter_", level, "_loess.png"), p, width=12, height=9, units="in")
      
      p <- ggplot(na.omit(dl[dl$trait_value>=0 & dl$range_value>=0,]),
                  aes(trait_value, range_value)) +
            geom_point(size=1.5, alpha=.3) +
            #geom_smooth(color="dodgerblue", se=F, size=1) +
            geom_smooth(color="red", se=F, size=1, method=lm) +
            facet_grid(range_name~trait_name, scales="free") +
            labs(title=paste0(title, ", raw\n"))
      ggsave(paste0("E:/BIEN/traits_ranges/v5/scatter_", level, ".png"), p, width=12, height=9, units="in")
      
      
      p <- ggplot(na.omit(dl[dl$trait_value>0 & dl$range_value>0,]),
                  aes(trait_value, range_value)) +
            geom_point(size=1.5, alpha=.3) +
            geom_smooth(color="dodgerblue", se=F, size=1) +
            geom_smooth(color="red", se=F, size=1, method=lm) +
            facet_grid(range_name~trait_name, scales="free") +
            scale_x_log10() +
            scale_y_log10() +
            labs(title=paste0(title, ", log-log\n"))
      #ggsave(paste0("E:/BIEN/traits_ranges/v5/scatter_loglog_", level, ".png"), p, width=12, height=9, units="in")
      
      dtl <- dl %>%
            group_by(trait_name) %>%
            mutate(percentile = ecdf(trait_value)(trait_value)) %>%
            filter(percentile > .025 & percentile < .975) %>%
            group_by(range_name) %>%
            mutate(percentile = ecdf(range_value)(range_value)) %>%
            filter(percentile > .025 & percentile < .975) %>%
            select(-percentile)
      
      p <- ggplot(dtl, aes(trait_value, range_value)) +
            geom_point(size=1.5, alpha=.3) +
            geom_smooth(color="dodgerblue", se=F, size=1) +
            geom_smooth(color="red", se=F, size=1, method=lm) +
            facet_grid(range_name~trait_name, scales="free") +
            scale_x_log10() +
            scale_y_log10() +
            labs(title=paste0(title, ", median 95%, log-log\n"))
      ggsave(paste0("E:/BIEN/traits_ranges/v5/scatter_trimmed_loglog_", level, ".png"), p, width=12, height=9, units="in")
      
      p <- ggplot(dtl[dtl$range_name=="updArea",], aes(trait_value, range_value)) +
            geom_point(size=1.5, alpha=.3) +
            geom_smooth(color="dodgerblue", se=F, size=1) +
            geom_smooth(color="red", se=F, size=1, method=lm) +
            facet_wrap(~trait_name, scales="free", ncol=4) +
            scale_x_log10() +
            scale_y_log10() +
            labs(title=paste0(title, ", median 95%, log-log, updArea\n"))
      #ggsave(paste0("E:/BIEN/traits_ranges/v5/scatter_trimmed_loglog_updArea_", level, ".png"), p, width=12, height=9, units="in")
      
      
      p <- ggplot(dtl, aes(trait_value, range_value)) +
            geom_point(size=1.5, alpha=.3) +
            geom_smooth(color="dodgerblue", se=F, size=1) +
            geom_smooth(color="red", se=F, size=1, method=lm) +
            facet_grid(range_name~trait_name, scales="free") +
            labs(title=paste0(title, ", median 95%\n"))
      #ggsave(paste0("E:/BIEN/traits_ranges/v5/scatter_trimmed_", level, ".png"), p, width=12, height=9, units="in")
      
}

######### same but by HABIT #########
level <- "species"

dl <- d
title <- paste0(level, "-level")

p <- ggplot(na.omit(dl[dl$trait_value>=0 & dl$range_value>=0,]),
            aes(trait_value, range_value, color=habit)) +
      geom_point(size=1.5, alpha=.3) +
      geom_smooth(se=F, size=1, method=lm) +
      facet_grid(range_name~trait_name, scales="free") +
      labs(title=paste0(title, ", raw\n")) +
      theme(legend.position="top")
ggsave(paste0("E:/BIEN/traits_ranges/v5/scatter_", level, "_HABIT.png"), p, width=12, height=9, units="in")

dtl <- dl %>%
      group_by(trait_name, habit) %>%
      mutate(percentile = ecdf(trait_value)(trait_value)) %>%
      filter(percentile > .025 & percentile < .975) %>%
      group_by(range_name, habit) %>%
      mutate(percentile = ecdf(range_value)(range_value)) %>%
      filter(percentile > .025 & percentile < .975) %>%
      select(-percentile)

p <- ggplot(dtl, aes(trait_value, range_value, color=habit)) +
      geom_point(size=1.5, alpha=.3) +
      geom_smooth(se=F, size=1, method=lm) +
      facet_grid(range_name~trait_name, scales="free") +
      scale_x_log10() +
      scale_y_log10() +
      labs(title=paste0(title, ", median 95%, log-log\n")) +
      theme(legend.position="top")
ggsave(paste0("E:/BIEN/traits_ranges/v5/scatter_trimmed_loglog_", level, "_HABIT.png"), p, width=12, height=9, units="in")

# range size by habit
dwl <- spread(dl, trait_name, trait_value)
p <- ggplot(dwl, aes(habit, range_value, color=habit, fill=habit)) +
      geom_boxplot(alpha=.5) +
      scale_y_log10() +
      facet_wrap(~range_name, scales="free", ncol=4) +
      labs(title="Range size by habit", y="range value (log scale)") +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            legend.position="top")
ggsave(paste0("E:/BIEN/traits_ranges/v5/range_by_habit_boxplot.png"), p, width=12, height=9, units="in")




############ within-family ###########

dl <- split(d, paste(d$family, d$trait_name, d$range_name))
scaleLog <- function(x){
      x$trait_value <- log10(x$trait_value)
      x$range_value <- log10(x$range_value)
      x$trait_value <- (x$trait_value - mean(x$trait_value, na.rm=T)) / sd(x$trait_value)
      x$range_value <- (x$range_value - mean(x$range_value, na.rm=T)) / sd(x$range_value)
      return(x)
}
dl <- do.call("rbind", lapply(dl, scaleLog))

level <- "within-family"
title <- "within-family anomalies, species-level"

p <- ggplot(dl,
            aes(trait_value, range_value)) +
      geom_point(size=1.5, alpha=.3) +
      geom_smooth(color="dodgerblue", se=F, size=1) +
      geom_smooth(color="red", se=F, size=1, method=lm) +
      facet_grid(range_name~trait_name, scales="free") +
      labs(title=paste0(title, ", rescaled log-log\n"))
ggsave(paste0("E:/BIEN/traits_ranges/v5/scatter_loglog_", level, ".png"), p, width=12, height=9, units="in")

p <- ggplot(dl[which((abs(dl$trait_value)<2) & (abs(dl$range_value)<2)),],
            aes(trait_value, range_value)) +
      geom_point(size=1.5, alpha=.3) +
      geom_smooth(color="dodgerblue", se=F, size=1) +
      geom_smooth(color="red", se=F, size=1, method=lm) +
      facet_grid(range_name~trait_name, scales="free") +
      labs(title=paste0(title, ", trimmed 2sd+, rescaled log-log\n"))
ggsave(paste0("E:/BIEN/traits_ranges/v5/scatter_trimmed_loglog_", level, ".png"), p, width=12, height=9, units="in")

### with HABIT ###
dl <- split(d, paste(d$family, d$trait_name, d$range_name, d$habit))
dl <- do.call("rbind", lapply(dl, scaleLog))

p <- ggplot(dl,
            aes(trait_value, range_value, color=habit)) +
      geom_point(size=1.5, alpha=.3) +
      geom_smooth(se=F, size=1, method=lm) +
      facet_grid(range_name~trait_name, scales="free") +
      labs(title=paste0(title, ", rescaled log-log\n"))
ggsave(paste0("E:/BIEN/traits_ranges/v5/scatter_loglog_", level, "_HABIT.png"), p, width=12, height=9, units="in")

p <- ggplot(dl[which((abs(dl$trait_value)<2) & (abs(dl$range_value)<2)),],
            aes(trait_value, range_value, color=habit)) +
      geom_point(size=1.5, alpha=.3) +
      geom_smooth(se=F, size=1, method=lm) +
      facet_grid(range_name~trait_name, scales="free") +
      labs(title=paste0(title, ", trimmed 2sd+, rescaled log-log\n"))
ggsave(paste0("E:/BIEN/traits_ranges/v5/scatter_trimmed_loglog_", level, "_HABIT.png"), p, width=12, height=9, units="in")


scaleRank <- function(x){
      x$trait_value <- scales::rescale(ecdf(x$trait_value)(x$trait_value))
      x$range_value <- scales::rescale(ecdf(x$range_value)(x$range_value))
      return(x)
}
dl <- split(d, paste(d$family, d$trait_name, d$range_name, d$habit))
dl <- do.call("rbind", lapply(dl, scaleRank))

p <- ggplot(dl,
            aes(trait_value, range_value, color=habit)) +
      geom_point(size=1.5, alpha=.3) +
      geom_smooth(se=F, size=1, method=lm) +
      facet_grid(range_name~trait_name, scales="free") +
      labs(title=paste0(title, ", rank-rescaled\n"), x="rank", y="rank")
ggsave(paste0("E:/BIEN/traits_ranges/v5/scatter_rank_", level, "_HABIT.png"), p, width=12, height=9, units="in")





############### experiment with regression ############

w <- d %>%
      mutate(trait_name=substr(trait_name, 1, 18)) %>%
      spread(trait_name, trait_value) %>%
      spread(range_name, range_value) %>%
      filter(habit=="Tree") %>%
      select(Leaf_Nmass, Specific_leaf_area, wood_density, updArea) %>%
      mutate_each(funs(log10), Leaf_Nmass, Specific_leaf_area) %>%
      mutate(updArea=sqrt(updArea)) %>%
      na.omit()

m <- lm(updArea ~ Leaf_Nmass + Specific_leaf_area + wood_density, data=w)
m <- lm(updArea ~ . + .*., data=w) # interactions

ww <- gather(w, variable, value, -updArea)
ggplot(ww, aes(value, updArea)) + geom_point() + geom_smooth(method=lm, se=F) +
      facet_wrap(~variable, scales="free")

dd <- split(d, paste(d$range_name, d$trait_name, sep="__"))
model <- function(x){
      require(caret)
      #x <- dd[[1]]
      x <- select(x, habit, trait_value, range_value) %>% 
            na.omit()
      trans <- function(y) predict(preProcess(y, method="BoxCox"), y)
      x <- cbind(x[,1], apply(x[,2:ncol(x)], 2, trans))
      
      m <- try(lm(range_value ~ . + .*., data=x))
      return(m)
}
dd <- lapply(dd, model)

### correlations

library(corrplot)
w <- d %>%
      mutate(trait_name=substr(trait_name, 1, 18)) %>%
      spread(trait_name, trait_value) %>%
      spread(range_name, range_value) %>%
      filter(habit=="Tree") %>%
      select(-family, -genus, -taxon, -habit) %>%
      log10() %>%
      as.matrix()
w[!is.finite(w)] <- NA
corrplot(cor(w, use="pairwise.complete.obs"))


#### random forest
library(randomForest)
library(missForest)

w <- d %>%
      mutate(trait_name=substr(trait_name, 1, 18)) %>%
      spread(trait_name, trait_value) %>%
      spread(range_name, range_value) %>%
      filter(habit=="Tree") %>%
      select(-family, -genus, -taxon, -habit, -sampleSize, -clippedHullArea)

mf <- missForest(as.matrix(w[,1:(ncol(w)-1)]))
rf <- randomForest(x=mf$ximp, y=w$updArea, importance=T)
varImpPlot(rf)




### next steps ###

# experiment with box-cox and rank rather than log transformations
# geographic/climatic subgroups
# regression to predict range size based on traits
# random forest for the same



