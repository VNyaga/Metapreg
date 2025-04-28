
#===========================
# sort dataframe by rmse
order_B<- unique(rbias[rbias$inf=="B" & rbias$Stat=="Median",] %>% arrange( desc(rel_bias)) %>% pull(id))
order_F <- unique(rbias[rbias$inf=="F",] %>% arrange(desc(rel_bias)) %>% pull(id))

# Plot
B <-  ggplot(data= simor[simor$inf=="B" & simor$Stat=="Median",]) +
  
  geom_half_violin(aes(x=factor(id, levels=order_B), y=esthat),
                   linewidth=0.2, side="r",  scale="area", color="gray", fill="gray") +
  
  geom_half_violin(data= simor[simor$inf=="B" & 
                                 simor$Covariance == orange  & 
                                 simor$Dist == "BN" &
                                 simor$Design == "comparative",], 
                   aes(x=factor(id, levels=order_B), y=esthat),
                   linewidth=0.2, side="r",  scale="area", color="orange", fill="orange") +
  
  geom_boxplot(aes(x=as.factor(id), y=esthat),
               width=0.1, color="black", alpha=0.2) +
  
  geom_text(data=coverage[coverage$inf=="B" & coverage$CImethod =='eti',], 
            aes(x=as.factor(id), y=right3, label=coverage), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$inf=="B" & coverage$CImethod =='hpd',], 
            aes(x=as.factor(id), y=right4, label=coverage), nudge_x = 0.2) +
  
  geom_text(data=rbias[rbias$inf=="B" & rbias$CImethod =='hpd' & rbias$Stat=="Median",], 
            aes(x=as.factor(id), y=right1, label=rel_bias), nudge_x = 0.2) +
  
  theme(legend.position="none") +
  
  coord_flip() +
  
  geom_hline(yintercept = trueor, linewidth=0.5, colour='black')   +

scale_y_continuous(position = 'right',
                   limits=c(minX, maxX),
                   breaks=c(minX, 
                            trueor,
                            min(rbias$right1),
                            min(coverage$right3),
                            min(coverage$right4),
                            maxX),
                   expand = c(0.06, 0.02),
                   labels = c("",  "True OR", "RBias", "Coverage(eti, hpd)", "", ""),
                   sec.axis = sec_axis(~.,breaks=c(trueor),
                                       labels=c(formatted_trueor))) +
  theme(axis.text.y=element_text(color='black'),
        axis.ticks.x=element_blank(),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'),
        axis.text.x=element_text(size=9),
        strip.text.x=element_text(size=9),
        panel.grid.major.y = element_line(color = 'gray'),
        panel.background = element_rect(fill = "white")) +
  xlab("") + ylab("") 

B
#=================================
ggsave('B-Sim.png',
       path = graphpath,
       width = 7,
       height = 3.5,
       dpi = 300)


#=========================================
# Plot
F <-  ggplot(data= simor[simor$inf=="F" & simor$Stat=="Mean",]) +
  
  geom_half_violin(aes(x=factor(id, levels=order_F), y=esthat),
                   linewidth=0.2, side="r",  scale="area", color="gray", fill="gray") +
  
  geom_half_violin(data= simor[simor$inf=="F" & 
                                 simor$Covariance == orange  & 
                                 simor$Dist == "BN" &
                                 simor$Design == "comparative" & 
                                 simor$Link=='LGT',], 
                   aes(x=factor(id, levels=order_F), y=esthat),
                   linewidth=0.2, side="r",  scale="area", color="orange", fill="orange") +
  
  geom_boxplot(aes(x=as.factor(id), y=esthat),
               width=0.1, color="black", alpha=0.2) +
  
  geom_text(data=coverage[coverage$inf=="F" & coverage$CImethod =='t',], 
            aes(x=as.factor(id), y=right3, label=coverage), nudge_x = 0.2) +
  
  geom_text(data=coverage[coverage$inf=="F" & coverage$CImethod =='z',], 
            aes(x=as.factor(id), y=right4, label=coverage), nudge_x = 0.2) +
  
  geom_text(data=rbias[rbias$inf=="F" & rbias$CImethod =='z',], 
            aes(x=as.factor(id), y=right1, label=rel_bias), nudge_x = 0.2) +
  
  theme(legend.position="none") +
  
  coord_flip() +
  
  facet_grid(cols=vars(Link)) +
  
  geom_hline(yintercept = trueor, linewidth=0.5, colour='black')   +
  
  scale_y_continuous(position = 'right',
                     limits=c(minX, maxX),
                     breaks=c(minX, 
                              trueor,
                              min(rbias$right1),
                              min(coverage$right3),
                              min(coverage$right4),
                              maxX),
                     expand = c(0.06, 0.02),
                     labels = c("",  "True OR", "RBias", "Coverage(t,z)", "", ""),
                     sec.axis = sec_axis(~.,breaks=c(trueor),
                                         labels=c(formatted_trueor))) +
  theme(axis.text.y=element_text(color='black'),
        axis.ticks.x=element_blank(),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'),
        strip.text.x=element_text(size=9),
        axis.text.x=element_text(size=9),
        panel.grid.major.y = element_line(color = 'gray'),
        panel.background = element_rect(fill = "white")) +
  xlab("") + ylab("")

F
#=================================
ggsave('F-Sim.png',
       path = graphpath,
       width = 10.5,
       height = 3.5,
       dpi = 300)

