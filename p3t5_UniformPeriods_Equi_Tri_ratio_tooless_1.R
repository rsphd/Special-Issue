library(MASS)
library(pracma)

#p: Number of periods, t: Number of treatments, n: Number of subjects
p=3
t=5
n=20

#Incidence matrix of period versus direct treatment effect for each subject
T_d1=matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0),nrow=p,ncol=t,byrow=T)
T_d2=matrix(c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0),nrow=p,ncol=t,byrow=T)
T_d3=matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1),nrow=p,ncol=t,byrow=T)
T_d4=matrix(c(0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0),nrow=p,ncol=t,byrow=T)
T_d5=matrix(c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0),nrow=p,ncol=t,byrow=T)
T_d6=matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0),nrow=p,ncol=t,byrow=T)
T_d7=matrix(c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0),nrow=p,ncol=t,byrow=T)
T_d8=matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1),nrow=p,ncol=t,byrow=T)
T_d9=matrix(c(0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0),nrow=p,ncol=t,byrow=T)
T_d10=matrix(c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0),nrow=p,ncol=t,byrow=T)
T_d11=matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0),nrow=p,ncol=t,byrow=T)
T_d12=matrix(c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0),nrow=p,ncol=t,byrow=T)
T_d13=matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1),nrow=p,ncol=t,byrow=T)
T_d14=matrix(c(0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0),nrow=p,ncol=t,byrow=T)
T_d15=matrix(c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0),nrow=p,ncol=t,byrow=T)
T_d16=matrix(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0),nrow=p,ncol=t,byrow=T)
T_d17=matrix(c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0),nrow=p,ncol=t,byrow=T)
T_d18=matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1),nrow=p,ncol=t,byrow=T)
T_d19=matrix(c(0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0),nrow=p,ncol=t,byrow=T)
T_d20=matrix(c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0),nrow=p,ncol=t,byrow=T)

#Combining period versus direct treatment effect incidence matrices for all subjects
T_d=rbind(T_d1, T_d2, T_d3, T_d4, T_d5, T_d6, T_d7, T_d8, T_d9, T_d10, T_d11, T_d12, T_d13, T_d14, T_d15, T_d16, T_d17, T_d18, T_d19, T_d20)

#Incidence matrix of period versus direct treatment effect for the first subject for orthogonal array design of Type I and strength 2 
T_d_ortho_1=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0),nrow=p,ncol=t,byrow=T)

psi=matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0),nrow=p,ncol=p,byrow=T)

#Identity matrices of dimension pxp, nxn and txt
I_p=diag(p)
I_n=diag(n)
I_t=diag(t)

#pxp, nxn and txt matrices with all elements as 1
J_p=matrix(rep(1,p^2),nrow=p,ncol=p,byrow=T)
J_n=matrix(rep(1,n^2),nrow=n,ncol=n,byrow=T)
J_t=matrix(rep(1,t^2),nrow=t,ncol=t,byrow=T)

H_n=I_n-J_n/n
H_t=I_t-J_t/t

#Incidence matrix of period versus carryover effect
F_d=kronecker(I_n,psi)%*%T_d

#Vector containing various values of r_(1)
r1=c(seq(-0.49,-0.01,0.01), seq(0.01,0.70,0.01))

Nr=c()
Dr=c()
efficiency=c()

#For loop to evaluate the efficiency of binary design at various values of r_(1) 
for(i in 1:length(r1))
{
  V_1=matrix(c(1,r1[i],r1[i],r1[i],1,r1[i],r1[i],r1[i],1),nrow=p,ncol=p,byrow=T)
  V_2=matrix(c(1,r1[i],0,r1[i],1,r1[i],0,r1[i],1),nrow=p,ncol=p,byrow=T)
  
  x=0.001 #sigma_1^2/sigma_2^2
  
  delta_1=1/sum(rowSums(solve(V_1)))
  V_1_star=solve(V_1)-delta_1*solve(V_1)%*%J_p%*%solve(V_1)
  A_1_star=kronecker(H_n,V_1_star)
  
  
  delta_2=1/sum(rowSums(solve(V_2)))
  V_2_star=solve(V_2)-delta_2*solve(V_2)%*%J_p%*%solve(V_2)
  A_2_star=kronecker(H_n,V_2_star)
  
  C_d11_1=t(T_d)%*%A_1_star%*%T_d
  C_d12_1=t(T_d)%*%A_1_star%*%F_d
  C_d21_1=t(C_d12_1)
  C_d22_1=t(F_d)%*%A_1_star%*%F_d
  
  C_d_1 = C_d11_1 - C_d12_1%*%pinv(C_d22_1)%*%C_d21_1
  
  C_d11_2=t(T_d)%*%A_2_star%*%T_d
  C_d12_2=t(T_d)%*%A_2_star%*%F_d
  C_d21_2=t(C_d12_2)
  C_d22_2=t(F_d)%*%A_2_star%*%F_d
  
  C_d_2 = C_d11_2 - C_d12_2%*%pinv(C_d22_2)%*%C_d21_2
  
  s1=sum(diag(C_d_1))
  s2=sum(diag(C_d_2))
  
  Nr[i]= s1+x*s2
  
  q_1_11 = sum(diag(t(T_d_ortho_1)%*%V_1_star%*%T_d_ortho_1))
  q_1_12 = sum(diag(t(T_d_ortho_1)%*%V_1_star%*%psi%*%T_d_ortho_1))
  q_1_22 = sum(diag(t(T_d_ortho_1)%*%t(psi)%*%V_1_star%*%psi%*%T_d_ortho_1)) - V_1_star[1,1]/t
  
  Q_1 = n/(t-1)*matrix(c(q_1_11,q_1_12,q_1_12,q_1_22),nrow=2,ncol=2,byrow=T)
  
  C_d_ortho_1 = (det(Q_1)/Q_1[2,2])*H_t
  
  q_2_11 = sum(diag(t(T_d_ortho_1)%*%V_2_star%*%T_d_ortho_1))
  q_2_12 = sum(diag(t(T_d_ortho_1)%*%V_2_star%*%psi%*%T_d_ortho_1))
  q_2_22 = sum(diag(t(T_d_ortho_1)%*%t(psi)%*%V_2_star%*%psi%*%T_d_ortho_1)) - V_2_star[1,1]/t
  
  Q_2 = n/(t-1)*matrix(c(q_2_11,q_2_12,q_2_12,q_2_22),nrow=2,ncol=2,byrow=T)
  
  C_d_ortho_2 = (det(Q_2)/Q_2[2,2])*H_t
  
  s1_ortho=sum(diag(C_d_ortho_1))
  s2_ortho=sum(diag(C_d_ortho_2))
  
  Dr[i]=s1_ortho+x*s2_ortho
  
  efficiency[i]=Nr[i]/Dr[i]
  
}

data=data.frame(efficiency)
names(data)=c("Efficiency")

#Calling library writexl
library(writexl)
#Storing dataframe in given path as excel file
write_xlsx(data,"D:/Special Issue/p3t5_UniformPeriods_Equi_Tri_ratio_tooless_1.xlsx")

#Calling library readxl
library(readxl)

p3t5_UniformPeriods_Equi_Tri_ratio_tooless_1=read_excel("D:/Special Issue/p3t5_UniformPeriods_Equi_Tri_ratio_tooless_1.xlsx")

min(p3t5_UniformPeriods_Equi_Tri_ratio_tooless_1$Efficiency) #Minimum efficiency

max(p3t5_UniformPeriods_Equi_Tri_ratio_tooless_1$Efficiency) #Maximum efficiency



