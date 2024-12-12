import { defineConfig } from "vitepress";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  lang: "en-US",
  title: "Deep mutational scanning of rabies G (Pasteur strain)",
  description:
    "Pseudovirus deep mutational scanning to measure how rabies G mutations affect cell entry and escape from a panel of monoclonal antibodies",
  base: "/RABV_Pasteur_G_DMS/",
  appearance: false,
  themeConfig: {
    // https://vitepress.dev/reference/default-theme-config
    nav: [
      { text: "Home", link: "/" },
      { text: "Appendix", link: "/appendix", target: "_self" },
    ],
    socialLinks: [{ icon: "github", link: "https://github.com/dms-vep/RABV_Pasteur_G_DMS" }],
    footer: {
      message: 'Study by Arjun Aditham, Caelan Radford, Caleb, Carr, & <a href="https://jbloomlab.org">Jesse Bloom</a>',
    },
  },
});