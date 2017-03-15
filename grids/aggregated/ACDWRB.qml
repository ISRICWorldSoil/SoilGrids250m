<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="2.14.1-Essen" minimumScale="0" maximumScale="1e+08" hasScaleBasedVisibilityFlag="0">
  <pipe>
    <rasterrenderer opacity="1" alphaBand="-1" classificationMax="4" classificationMinMaxOrigin="User" band="1" classificationMin="0" type="singlebandpseudocolor">
      <rasterTransparency/>
      <rastershader>
        <colorrampshader colorRampType="INTERPOLATED" clip="0">
          <item alpha="255" value="0" label="Low or none acidity problems" color="#fff5eb"/>
          <item alpha="255" value="1" label="Slightly acid (pH&lt;6.6) and/or dystric soil properties" color="#fdd1a5"/>
          <item alpha="255" value="2" label="Moderately acid (pH&lt;5.5) and/or dystric Podzols and/or Gleysols" color="#fd9243"/>
          <item alpha="255" value="3" label="Strongly acid (pH&lt;5) and/or dystric Acrisols and/or Alisols" color="#de4f05"/>
          <item alpha="255" value="4" label="Extremely acid (pH&lt;4.5) and/or Acrisols, Alisols and/or alumic properties" color="#7f2704"/>
        </colorrampshader>
      </rastershader>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0"/>
    <huesaturation colorizeGreen="128" colorizeOn="0" colorizeRed="255" colorizeBlue="128" grayscaleMode="0" saturation="0" colorizeStrength="100"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
