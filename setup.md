---
layout: page
title: Setup
---

{% comment %} Setup {% endcomment %}
<h2 id="setup">Setup</h2>
<p>
  To participate in a this workshop, you will need access to UF Hipergator and the software described below.
</p>
<p>
  We maintain a list of common issues that occur during installation as a reference for instructors
  that may be useful on the
  <a href = "{{ site.baseSite }}{{ site.faq }}">FAQ page</a>.
</p>

{% comment %} Zoom {% endcomment %}
{% include setup/conference.html %}

{% comment %} UF VPN Access {% endcomment %}
<div id="vpn" markdown="1">
### Connecting to Gatorlink VPN
Gatorlink VPN is required to access University if Florida Research Computing OnDemand
service, which will provide access to HiperGator cluster.
1. Visit [UF VPN webpage (vpn.ufl.edu)](https://vpn.ufl.edu){: target="_blank"}.
2. Login using the username and password that was emailed to you.
3. In the next page, download **Cisco Anyconnect** software and install it.
4. Once installed, run **Cisco Anyconnect** and click **Connect**.
5. Enter the username and password provided and click **OK**. 
You should now be connected to Gatorlink VPN.
</div>

{% comment %} Hipergator Access {% endcomment %}
<div id="shell" markdown="1">
### Connecting to HiperGator
Make sure you are connected to Gatorlink VPN before the following steps.
1. In you web browser, navigate to [UFRC OnDemand (ood.rc.ufl.edu)](https://ood.rc.ufl.edu/){: target="_blank"}.
2. Login using provided username and password. You will be redirected to UFRC OnDemand homepage.
3. To connect to remote shell, click on **Clusters** in navbar and click 
**Hipergator Shell Access**. Note that the shell will open in new tab.
4. To view and transfer your files, click on **Files** in navbar 
in UFRC OnDemand homepage and click **/blue/general_workshop**.
Double click on directory named as your username to find your files.
</div>

{% comment %} Text editor {% endcomment %}
{% include setup/editor.html %}

{% comment %} Figtree {% endcomment %}
{% include setup/figtree.md %}