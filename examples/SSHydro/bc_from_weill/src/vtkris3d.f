         subroutine vtkris3d(ntetra, nnode,iunit,tetra,
     1                       u,v,vx,vy,vz,x,y,z,ks,time,vtkf)
         implicit none
         integer nnode, ntetra, iunit,  i, j, numb
         integer tetra(5,*), field,vtkf
         real*8 ks(*)
         real*8 u(*), v(*),x(*), y(*), z(*)
         real*8  vx(*),vy(*),vz(*)
         real*8 time
         character*15 nome
         character*15 file_vtk

         field=1

         write(nome,'(i3,a4)') iunit, '.vtk'
         file_vtk='vtk'//'/'//nome
         open(iunit, file=file_vtk)
cm       open(iunit, file=nome)
        
         write(iunit,78) 
 78      FORMAT('# vtk DataFile Version 2.0',
     1     /,'3D Unstructured Grid of Linear Triangles', 
     2     /,'ASCII')

         write(iunit,77) field, field, field, time
 77      FORMAT('DATASET UNSTRUCTURED_GRID',
     1     /, 'FIELD FieldData',1x, i2,
     2     /, 'TIME', 1x, i1,1x,i1, 1x, 'double',
     3     /, 1f18.5)

        write(iunit,79) nnode
 79     FORMAT('POINTS',1x, i8, ' float')

         do i=1,nnode
             write(iunit,100) x(i), y(i), z(i)
         end do
 100     format(4x, 3(1pe16.8))

         numb=ntetra*5
         write(iunit, 80) ntetra,numb
 80      format('CELLS',1x,i8,1x,i8)
         j=4
         do i=1,ntetra
           write(iunit,101) j, tetra(1,i)-1, tetra(2,i)-1, 
     1                    tetra(3,i)-1, tetra(4,i)-1
         end do
 101     format(i1, 3x,i8,3x,i8,3x,i8,3x,i8)

         write(iunit,81)  ntetra
 81      format('CELL_TYPES',i8)
         do i=1,ntetra
          write(iunit,*) 10
         end do

         write(iunit,82)  nnode
 82      format('POINT_DATA', 1x,i8, 
     1    /, 'SCALARS pressure float',
     2    /, 'LOOKUP_TABLE default' )
         do i=1,nnode
          write(iunit,*) u(i)
         end do
         if (vtkf .ge.2) then
         write(iunit,83) 
 83      format('SCALARS saturation float'
     1    /, 'LOOKUP_TABLE default' )
         do i=1,nnode
          write(iunit,*) v(i)
         end do
         end if

         if (vtkf.ge.3) then
         write(iunit,84)  ntetra
 84      format('CELL_DATA', 1x,i8, 
     1   /,'SCALARS permeability float',
     1    /, 'LOOKUP_TABLE default' )
         do i=1,ntetra
             write(iunit,*) ks(i)
         end do
         end if

         if (vtkf.ge.4) then
         write(iunit,85) 
 85      format( 'VECTORS velocity float')
         do i=1,ntetra
             write(iunit,*) vx(i), vy(i), vz(i)
         end do
         end if
         close(iunit)
         return
         end
         
